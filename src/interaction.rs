use itertools::Itertools;
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::prelude::*;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Interaction {
    pub n: usize,
    consv: Vec<Vec<i64>>,
    pub edges: Vec<Vec<(usize, usize)>>,
}

impl Interaction {
    /// Construct the trivial (empty) interaction of size n
    fn new(n: usize) -> Self {
        let mut consv: Vec<Vec<i64>> = vec![vec![0; n]; n - 1];
        (1..n).for_each(|i| consv[i - 1][i] = 1);

        Self {
            n,
            consv,
            edges: vec![],
        }
    }

    /// Construct an interaction by given consv and edges.
    fn create_unchecked(n: usize, consv: Vec<Vec<i64>>, edges: Vec<Vec<(usize, usize)>>) -> Self {
        Self { n, consv, edges }
    }

    /// Get the list of edges of an interaction with the given conserved quantity.
    fn get_edges_from_consv(n: usize, consv: Vec<Vec<i64>>) -> Vec<Vec<(usize, usize)>> {
        let mut hm: HashMap<Vec<i64>, Vec<(usize, usize)>> = HashMap::new();

        // We classify verticies by values of conserved quantities.
        (0..n).combinations_with_replacement(2).for_each(|v| {
            let consv_vector: Vec<i64> = consv.iter().map(|xi| xi[v[0]] + xi[v[1]]).collect();

            match hm.get_mut(&consv_vector) {
                Some(x) => {
                    x.push((v[0], v[1]));
                }
                None => {
                    hm.insert(consv_vector.clone(), vec![(v[0], v[1])]);
                }
            }
        });

        let mut edges: Vec<Vec<(usize, usize)>> = vec![];

        // We add an edge if verticies have same values of conserved quantities.
        hm.into_iter()
            .filter(|(_, val)| val.len() > 1)
            .for_each(|(_, mut val)| {
                val.sort();
                edges.push(val);
            });

        // We always assume that this list is sorted.
        edges.sort();
        edges
    }

    /// Construct an interaction by given consv and edges.
    fn create_from_consv(n: usize, consv: Vec<Vec<i64>>) -> Self {
        Interaction::create_unchecked(
            n,
            consv.clone(),
            Interaction::get_edges_from_consv(n, consv),
        )
    }

    /// Check whether this interaction is separable.
    // In this function, we find a pair of states such that
    // these values of each conserved quantity are same.
    fn is_separable(&self) -> bool {
        let mut hs: HashSet<Vec<i64>> = HashSet::new();

        (0..self.n).all(|i| {
            let consv_values: Vec<i64> = self.consv.iter().map(|v| v[i]).collect();
            hs.insert(consv_values)
        })
    }

    /// Add an edge to the interaction
    fn merge(&self, (a, b): (usize, usize), (c, d): (usize, usize)) -> Option<Interaction> {
        let mut new_consv = vec![];

        // The algorihm is given in our paper [4, Lemma 4.1].
        let base_xi = self
            .consv
            .iter()
            .find(|xi| xi[a] + xi[b] != xi[c] + xi[d])?;

        let diff_b = base_xi[c] + base_xi[d] - base_xi[a] - base_xi[b];

        self.consv.iter().for_each(|xi| {
            let diff = xi[c] + xi[d] - xi[a] - xi[b];
            if xi != base_xi {
                let new_xi: Vec<i64> = xi
                    .iter()
                    .enumerate()
                    .map(|(i, &x)| diff_b * x - diff * base_xi[i])
                    .collect();
                new_consv.push(new_xi);
            }
        });

        Some(Interaction::create_from_consv(self.n, new_consv))
    }
}

pub struct InteractionsModEquiv {
    n: usize,
    // This is a hash set which saving the edge_list of interactions.
    my_inter_hs: HashSet<Vec<Vec<(usize, usize)>>>,
    // This is a list of interactions.
    my_inter_list: Vec<Interaction>,
    // This is a list of candidates of an edge of an interaction.
    new_edge_list: Vec<(usize, usize, usize, usize)>,
}

impl InteractionsModEquiv {
    /// Initialize our list.
    pub fn new(n: usize) -> Self {
        let mut new_edge_list: Vec<(usize, usize, usize, usize)> = vec![];
        (0..n).combinations_with_replacement(2).for_each(|origin| {
            (origin[0] + 1..n)
                .combinations_with_replacement(2)
                .filter(|target| origin[1] != target[0] && origin[1] != target[1])
                .for_each(|target| {
                    new_edge_list.push((origin[0], origin[1], target[0], target[1]));
                });
        });

        Self {
            n,
            my_inter_hs: HashSet::new(),
            my_inter_list: vec![],
            new_edge_list,
        }
    }

    /// Save the list of interactions.
    pub fn output_json(&self, file_name: String) -> std::io::Result<()> {
        let serialized: String = serde_json::to_string(&(self.my_inter_list)).unwrap();
        let mut file = File::create(file_name)?;
        file.write_all(serialized.as_bytes())?;
        Ok(())
    }

    /// Add an interaction to our list.
    fn add(&mut self, new_inter: Interaction) -> bool {
        let edges = new_inter.edges.clone();
        let consv = new_inter.consv.clone();

        // Check whether the interaction is weakly equivarent to the other interaction.
        if (0..self.n).permutations(self.n).all(|perm| {
            let perm_consv: Vec<Vec<i64>> = consv
                .iter()
                .map(|xi| perm.iter().map(|&x| xi[x]).collect())
                .collect();
            let permuted_interaction = Interaction::create_from_consv(self.n, perm_consv);
            !self.my_inter_hs.contains(&permuted_interaction.edges)
        }) {
            self.my_inter_hs.insert(edges);
            self.my_inter_list.push(new_inter);
            true
        } else {
            false
        }
    }

    /// The starter of our construction program.
    pub fn create_list(&mut self) {
        let trivial_inter = Interaction::new(self.n);
        self.add(trivial_inter.clone());
        self.add_to_list(trivial_inter, 0);
    }

    // Create list recursively......
    fn add_to_list(&mut self, inter: Interaction, index: usize) {
        if inter.consv.len() == 1 {
            return;
        }
        for i in index..self.new_edge_list.len() {
            let (a, b, c, d) = self.new_edge_list[i];
            let Some(new_inter) = inter.merge((a, b), (c, d)) else {
                continue;
            };
            if new_inter.is_separable() && self.add(new_inter.clone()) {
                self.add_to_list(new_inter, i + 1);
            }
        }
    }
}
