//! Coupling-owned spatial records for distributed SPH-CFD exchange.

use field_core::PartitionDirectory;
use grass_multi::{EntityId, ReceivedPayload, RoutedPayload};
use std::collections::BTreeMap;

const PARCEL_BYTES: usize = 8 + 4 + 3 * 8 + 3 * 8 + 8;
const FORCE_BYTES: usize = 8 + 4 + 3 * 8;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RoutedParcel {
    pub id: u64,
    pub sph_owner: i32,
    pub center: [f64; 3],
    pub velocity: [f64; 3],
    pub volume: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RoutedParcelForce {
    pub id: u64,
    pub sph_owner: i32,
    pub force: [f64; 3],
}

pub fn route_parcels(
    directory: &PartitionDirectory,
    parcels: &[RoutedParcel],
) -> Result<Vec<RoutedPayload>, String> {
    parcels
        .iter()
        .map(|parcel| {
            let destination = directory
                .owner_rank(parcel.center)
                .ok_or_else(|| format!("parcel {} is outside the CFD domain", parcel.id))?;
            Ok(RoutedPayload::new(
                destination,
                EntityId(parcel.id),
                encode_parcel(*parcel),
            ))
        })
        .collect()
}

pub fn decode_parcels(records: &[ReceivedPayload]) -> Result<Vec<RoutedParcel>, String> {
    records
        .iter()
        .map(|record| {
            let parcel = decode_parcel(&record.payload)?;
            (parcel.id == record.entity_id.0)
                .then_some(parcel)
                .ok_or_else(|| "parcel entity ID does not match its frame".to_owned())
        })
        .collect()
}

pub fn route_forces(forces: &[RoutedParcelForce]) -> Vec<RoutedPayload> {
    forces
        .iter()
        .map(|force| RoutedPayload::new(force.sph_owner, EntityId(force.id), encode_force(*force)))
        .collect()
}

pub fn reduce_forces(records: &[ReceivedPayload]) -> Result<Vec<RoutedParcelForce>, String> {
    let mut reduced = BTreeMap::<u64, RoutedParcelForce>::new();
    for record in records {
        let force = decode_force(&record.payload)?;
        if force.id != record.entity_id.0 {
            return Err("force entity ID does not match its frame".to_owned());
        }
        reduced
            .entry(force.id)
            .and_modify(|sum| {
                for axis in 0..3 {
                    sum.force[axis] += force.force[axis];
                }
            })
            .or_insert(force);
    }
    Ok(reduced.into_values().collect())
}

fn encode_parcel(parcel: RoutedParcel) -> Vec<u8> {
    let mut bytes = Vec::with_capacity(PARCEL_BYTES);
    bytes.extend_from_slice(&parcel.id.to_le_bytes());
    bytes.extend_from_slice(&parcel.sph_owner.to_le_bytes());
    for value in parcel.center.into_iter().chain(parcel.velocity) {
        bytes.extend_from_slice(&value.to_le_bytes());
    }
    bytes.extend_from_slice(&parcel.volume.to_le_bytes());
    bytes
}

fn decode_parcel(bytes: &[u8]) -> Result<RoutedParcel, String> {
    if bytes.len() != PARCEL_BYTES {
        return Err(format!(
            "parcel payload has {} bytes, expected {PARCEL_BYTES}",
            bytes.len()
        ));
    }
    let mut at = 0;
    let id = take_u64(bytes, &mut at);
    let sph_owner = take_i32(bytes, &mut at);
    let center = std::array::from_fn(|_| take_f64(bytes, &mut at));
    let velocity = std::array::from_fn(|_| take_f64(bytes, &mut at));
    let volume = take_f64(bytes, &mut at);
    Ok(RoutedParcel {
        id,
        sph_owner,
        center,
        velocity,
        volume,
    })
}

fn encode_force(force: RoutedParcelForce) -> Vec<u8> {
    let mut bytes = Vec::with_capacity(FORCE_BYTES);
    bytes.extend_from_slice(&force.id.to_le_bytes());
    bytes.extend_from_slice(&force.sph_owner.to_le_bytes());
    for value in force.force {
        bytes.extend_from_slice(&value.to_le_bytes());
    }
    bytes
}

fn decode_force(bytes: &[u8]) -> Result<RoutedParcelForce, String> {
    if bytes.len() != FORCE_BYTES {
        return Err("invalid force payload size".to_owned());
    }
    let mut at = 0;
    Ok(RoutedParcelForce {
        id: take_u64(bytes, &mut at),
        sph_owner: take_i32(bytes, &mut at),
        force: std::array::from_fn(|_| take_f64(bytes, &mut at)),
    })
}

fn take_u64(bytes: &[u8], at: &mut usize) -> u64 {
    let value = u64::from_le_bytes(bytes[*at..*at + 8].try_into().unwrap());
    *at += 8;
    value
}
fn take_i32(bytes: &[u8], at: &mut usize) -> i32 {
    let value = i32::from_le_bytes(bytes[*at..*at + 4].try_into().unwrap());
    *at += 4;
    value
}
fn take_f64(bytes: &[u8], at: &mut usize) -> f64 {
    let value = f64::from_le_bytes(bytes[*at..*at + 8].try_into().unwrap());
    *at += 8;
    value
}
