use crate::create_overlap_graph::OverlapGraph;
use flate2::read::MultiGzDecoder;
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::Path;

/// Opens a FASTQ file for buffered reading, automatically decompressing if gzipped.
pub fn open_fastq_reader(path: &Path) -> std::io::Result<Box<dyn BufRead>> {
    let mut file = File::open(path)?;
    let mut magic = [0u8; 2];
    let is_gz = match file.read_exact(&mut magic) {
        Ok(()) => magic == [0x1f, 0x8b],
        Err(_) => false,
    };
    file.seek(SeekFrom::Start(0))?;

    if is_gz {
        let gz = MultiGzDecoder::new(BufReader::with_capacity(128 * 1024, file));
        Ok(Box::new(BufReader::with_capacity(128 * 1024, gz)))
    } else {
        Ok(Box::new(BufReader::with_capacity(128 * 1024, file)))
    }
}

/// Get the reverse-complement of a node (flip trailing '+' <-> '-').
pub fn rc_node(id: &str) -> String {
    if let Some(last) = id.chars().last() {
        if last == '+' {
            let base = &id[..id.len() - 1];
            return format!("{}-", base);
        } else if last == '-' {
            let base = &id[..id.len() - 1];
            return format!("{}+", base);
        }
    }
    id.to_string()
}

/// Delete a set of nodes (both orientations) from the graph and remove associated edges
pub fn delete_nodes_and_edges(graph: &mut OverlapGraph, nodes_to_delete: &HashSet<String>) {
    // Initialize set of nodes to remove
    let mut oriented_nodes_to_delete: HashSet<String> = HashSet::new();

    for node in nodes_to_delete.iter() {
        // add to set of nodes to delete
        oriented_nodes_to_delete.insert(node.clone());
        // add rc counterpart
        oriented_nodes_to_delete.insert(rc_node(node));
    }

    // Delete nodes from graph.nodes
    for oriented_node in oriented_nodes_to_delete.iter() {
        graph.nodes.remove(oriented_node);
    }

    // Remove edges pointing to removed nodes
    for (_src, node) in graph.nodes.iter_mut() {
        node.edges
            .retain(|e| !oriented_nodes_to_delete.contains(&e.target_id));
    }
}

pub fn rev_comp(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' | 'a' => 'T',
            'T' | 't' => 'A',
            'C' | 'c' => 'G',
            'G' | 'g' => 'C',
            'N' | 'n' => 'N',
            _ => c,
        })
        .collect()
}
