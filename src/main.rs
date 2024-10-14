use itertools::Itertools;
use proconio::input;
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::prelude::*;

#[derive(Clone, Debug, Serialize, Deserialize)]
struct Interaction {
    n: usize,
    consv: Vec<Vec<i64>>,
    edges: Vec<Vec<(usize, usize)>>,
}

impl Interaction {
    // construct the trivial (empty) interaction of size n
    fn new(n: usize) -> Self {
        let mut consv: Vec<Vec<i64>> = vec![vec![0; n]; n - 1];
        (1..n).for_each(|i| consv[i - 1][i] = 1);

        Self {
            n,
            consv,
            edges: vec![],
        }
    }

    // Check whether this interaction is separable.
    // In this function, we find a pair of states such that
    // these values of each conserved quantity are same.
    fn is_separable(&mut self, consv_list: Vec<Vec<i64>>) -> bool {
        // This is a hash set saving values of each conserved quantity.
        let mut hs: HashSet<Vec<i64>> = HashSet::new();

        // Loop over the set of states.
        (0..self.n).all(|i| {
            // This is a vector of values of conserved quantities of the state "i".
            let mut consv_values: Vec<i64> = vec![];
            consv_list.iter().for_each(|v| consv_values.push(v[i]));
            // for v in consv_list.iter() {
            // consv_values.push(v[i]);
            // }

            // Is there other state which has same values of conserved quantities?
            hs.insert(consv_values)
        })
    }

    // Add an edge to the interaction
    fn merge(&mut self, (a, b): (usize, usize), (c, d): (usize, usize)) -> bool {
        // If the dimension of the space of conserved quantities is 1,
        // we cannot add an edges not to trivialize the interaction.
        if self.consv.len() == 1 {
            return false;
        }

        let mut new_consv = vec![];

        // The algorihm is given in our paper.
        // let mut base_xi = vec![0; self.n];

        let Some(base_xi) = self.consv.iter().find(|xi| xi[a] + xi[b] != xi[c] + xi[d]) else {
            return false;
        };
        // for xi in self.consv.iter() {
        // if xi[a] + xi[b] != xi[c] + xi[d] {
        // base_xi = xi.clone();
        // flag = false;
        // }
        // }

        // If No, this edge are already contained in the interaction.
        // if flag {
        // return false;
        // }

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

        // for xi in self.consv.iter() {
        // let diff = xi[c] + xi[d] - xi[a] - xi[b];
        // if *xi != base_xi {
        // let new_xi: Vec<i64> = xi
        // .iter()
        // .enumerate()
        // .map(|(i, &x)| diff_b * x - diff * base_xi[i])
        // .collect();
        // new_consv.push(new_xi);
        // }
        // }

        // Is the new interaction separable?
        if !self.is_separable(new_consv.clone()) {
            false
        } else {
            self.consv = new_consv.clone();
            true
        }
    }

    // Get the list of edges of the interaction.
    fn get_edges(&mut self) -> Vec<Vec<(usize, usize)>> {
        let mut hm: HashMap<Vec<i64>, Vec<(usize, usize)>> = HashMap::new();

        // We classify verticies by values of conserved quantities.
        for v in (0..self.n).combinations_with_replacement(2) {
            let mut consv_vector: Vec<i64> = vec![];
            for xi in &self.consv {
                consv_vector.push(xi[v[0]] + xi[v[1]]);
            }

            match hm.get_mut(&consv_vector) {
                Some(x) => {
                    x.push((v[0], v[1]));
                }
                None => {
                    hm.insert(consv_vector.clone(), vec![(v[0], v[1])]);
                }
            }
        }

        let mut edges: Vec<Vec<(usize, usize)>> = vec![];

        // We add an edge if verticies have same values of conserved quantities.
        for (_, x) in hm.iter() {
            if x.len() > 1 {
                let mut x_sort = x.clone();
                x_sort.sort();
                edges.push(x_sort.clone());
            }
        }

        // We always assume that this list is sorted.
        edges.sort();
        self.edges = edges.clone();

        edges
    }
}

struct InteractionsModEquiv {
    n: usize,
    // This is a hash set which saving the edge_list of interactions.
    my_inter_hs: HashSet<Vec<Vec<(usize, usize)>>>,
    // This is a list of interactions.
    my_inter_list: Vec<Interaction>,
    // This is a list of candidates of an edge of an interaction.
    new_edge_list: Vec<(usize, usize, usize, usize)>,
}

impl InteractionsModEquiv {
    // Initialize our list.
    fn new(n: usize) -> Self {
        let mut new_edge_list: Vec<(usize, usize, usize, usize)> = vec![];
        for origin in (0..n).combinations_with_replacement(2) {
            let min = origin[0];
            for target in (min + 1..n).combinations_with_replacement(2) {
                if origin[1] != target[0] && origin[1] != target[1] {
                    new_edge_list.push((origin[0], origin[1], target[0], target[1]));
                }
            }
        }

        Self {
            n,
            my_inter_hs: HashSet::new(),
            my_inter_list: vec![],
            new_edge_list,
        }
    }

    // Save the list of interactions.
    fn output_json(&mut self, file_name: String) -> std::io::Result<()> {
        let serialized: String = serde_json::to_string(&(self.my_inter_list)).unwrap();
        let mut file = File::create(file_name)?;
        file.write_all(serialized.as_bytes())?;
        Ok(())
    }

    // Get the list of edges from the conserved quantities.
    // The algorithm is same as the edges_list function above.
    fn edges_list(&mut self, consv: Vec<Vec<i64>>) -> Vec<Vec<(usize, usize)>> {
        let mut hm: HashMap<Vec<i64>, Vec<(usize, usize)>> = HashMap::new();
        for v in (0..self.n).combinations_with_replacement(2) {
            let mut temp_vector: Vec<i64> = vec![];
            for xi in consv.iter() {
                temp_vector.push(xi[v[0]] + xi[v[1]]);
            }

            match hm.get_mut(&temp_vector) {
                Some(x) => {
                    x.push((v[0], v[1]));
                }
                None => {
                    hm.insert(temp_vector.clone(), vec![(v[0], v[1])]);
                }
            }
        }

        let mut edges: Vec<Vec<(usize, usize)>> = vec![];

        for (_, x) in hm.iter() {
            if x.len() > 1 {
                let mut x_sort = x.clone();
                x_sort.sort();
                edges.push(x_sort.clone());
            }
        }

        edges.sort();
        edges
    }

    // Add an interaction to our list.
    fn add(&mut self, mut new_inter: Interaction) -> bool {
        let edges = new_inter.get_edges();
        let consv = new_inter.consv.clone();

        // Check whether the interaction is weakly equivarent to the other interaction.
        for perm in (0..self.n).permutations(self.n) {
            let mut perm_consv: Vec<Vec<i64>> = vec![];
            for xi in consv.iter() {
                let mut perm_xi: Vec<i64> = vec![];
                for &x in perm.iter() {
                    perm_xi.push(xi[x]);
                }
                perm_consv.push(perm_xi.clone());
            }

            let edges_permuted = self.edges_list(perm_consv);

            match self.my_inter_hs.contains(&edges_permuted) {
                true => {
                    return false;
                }
                false => {}
            }
        }

        self.my_inter_list.push(new_inter);
        self.my_inter_hs.insert(edges.clone());
        true
    }

    // This is the starter of our construction program.
    fn create_list(&mut self) {
        let trivial_inter = Interaction::new(self.n);
        self.add(trivial_inter.clone());
        self.add_to_list(trivial_inter, 0);
    }

    // Create list recursively......
    fn add_to_list(&mut self, inter: Interaction, index: usize) {
        for i in index..self.new_edge_list.len() {
            let (a, b, c, d) = self.new_edge_list[i];
            let mut new_inter = inter.clone();
            if new_inter.merge((a, b), (c, d)) && self.add(new_inter.clone()) {
                self.add_to_list(new_inter.clone(), i + 1);
            }
        }
    }
}

fn main() {
    input! {
        n: usize,
    }

    let mut inter_list = InteractionsModEquiv::new(n);
    inter_list.create_list();
    let file_name = String::from("output/size_4.json");
    let result = inter_list.output_json(file_name);
    match result {
        Ok(..) => {}
        Err(err) => {
            println!("{err}");
        }
    }
}
