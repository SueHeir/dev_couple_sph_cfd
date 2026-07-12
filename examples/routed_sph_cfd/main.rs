//! The same binary and TOML run a local SPH-CFD route or a 3+2 MPI split.

use field_core::{PartitionDirectory, UniformMeshConfig};
use grass_multi::{CoupledPairRunner, CouplingEpoch, PairRun, RoleLaunch};
use sph_cfd::routing::{
    decode_parcels, reduce_forces, route_forces, route_parcels, RoutedParcel, RoutedParcelForce,
};

const DEFAULT_CONFIG: &str = include_str!("config.toml");

fn directory(parts: i32) -> PartitionDirectory {
    PartitionDirectory::from_uniform_config(
        &UniformMeshConfig {
            nx: 12,
            ny: 2,
            nz: 2,
            ng: 1,
            bounds_lo: [0.0; 3],
            bounds_hi: [1.0; 3],
            y_edges: None,
            z_edges: None,
        },
        [parts, 1, 1],
    )
}

fn load(parcel: RoutedParcel) -> [f64; 3] {
    let fluid_velocity = [parcel.center[0], 0.0, 0.25];
    std::array::from_fn(|axis| 1.5 * parcel.volume * (fluid_velocity[axis] - parcel.velocity[axis]))
}

fn run_role(launch: RoleLaunch) -> usize {
    let role = launch.role().to_owned();
    let exchange = launch.into_routed_exchange();
    let (rank, size) = exchange.role_position();
    let parcels = if role == "sph" {
        vec![RoutedParcel {
            id: 500 + rank as u64,
            sph_owner: rank,
            center: [(rank as f64 + 0.5) / size as f64, 0.5, 0.5],
            velocity: [0.02 * (rank + 1) as f64, 0.0, 0.0],
            volume: 1.0e-6,
        }]
    } else {
        Vec::new()
    };
    let outgoing = if role == "sph" {
        route_parcels(&directory(exchange.peer_size()), &parcels).expect("route SPH parcels")
    } else {
        Vec::new()
    };
    let incoming = exchange
        .exchange(CouplingEpoch(0), &outgoing)
        .expect("exchange parcels");
    let force_routes = if role == "cfd" {
        let owners = directory(size);
        let forces = decode_parcels(&incoming)
            .expect("decode parcels")
            .into_iter()
            .map(|parcel| {
                assert_eq!(owners.owner_rank(parcel.center), Some(rank));
                RoutedParcelForce {
                    id: parcel.id,
                    sph_owner: parcel.sph_owner,
                    force: load(parcel),
                }
            })
            .collect::<Vec<_>>();
        route_forces(&forces)
    } else {
        Vec::new()
    };
    let returned = exchange
        .exchange(CouplingEpoch(1), &force_routes)
        .expect("return SPH forces");
    if role == "sph" {
        let forces = reduce_forces(&returned).expect("reduce SPH forces");
        assert_eq!(forces.len(), parcels.len());
        assert_eq!(forces[0].force, load(parcels[0]));
        forces.len()
    } else {
        incoming.len()
    }
}

fn main() {
    let run = CoupledPairRunner::from_cli_or(DEFAULT_CONFIG)
        .and_then(|runner| runner.run(run_role, run_role))
        .unwrap_or_else(|error| panic!("run routed SPH-CFD example: {error}"));
    match run {
        PairRun::Local { first, second } => {
            println!("LOCAL sph_records={first} cfd_records={second}");
            println!("PASS same-binary local routed SPH-CFD coupling");
        }
        PairRun::Split { role, result } => println!("MPI role={role} local_records={result}"),
    }
}
