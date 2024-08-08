use proconio::input;
use std::collections::{HashMap, HashSet};
use itertools::Itertools;
use serde::{Serialize, Deserialize};
use std::fs::File;
use std::io::prelude::*;

#[derive(Clone, Debug)]
#[derive(Serialize, Deserialize)]
struct Interaction {
    n: usize,
    consv: Vec<Vec<i64>>,
    edges: Vec<Vec<(usize,usize)>>,
}

impl Interaction {

    fn new(n: usize) -> Self {
        let mut consv: Vec<Vec<i64>> = vec![vec![0;n];n-1];
        for i in 1..n {
            consv[i-1][i] = 1;
        }

        Self {
            n, 
            consv,
            edges: vec![],
        }
    }

    fn is_separable(&mut self, consv_list: Vec<Vec<i64>>) -> bool {

        let mut hs: HashSet<Vec<i64>> = HashSet::new();

        for i in 0..self.n {
            let mut temp_vector:Vec<i64> = vec![];
            for j in 0..consv_list.len() {
                temp_vector.push(consv_list[j][i]);
            }

            match hs.contains(&temp_vector) {
                true => {
                    return true
                }, 
                false => {
                    hs.insert(temp_vector.clone());
                }
            }
        }
        return false;
    }

    fn merge(&mut self, (a,b): (usize,usize), (c,d): (usize,usize)) -> bool {
        
        if self.consv.len() == 1 {
            return false; 
        }
        let mut new_consv = vec![];
        let mut base_xi = vec![0;self.n];

        let mut flag = true;
        for xi in self.consv.iter() {
            if xi[a] + xi[b] != xi[c] + xi[d] {
                base_xi = xi.clone();
                flag = false;
            }
        }
        
        if flag {
            return false;
        }

        let diff_b = base_xi[c] + base_xi[d] - base_xi[a] - base_xi[b];

        for xi in self.consv.iter() {
            let diff = xi[c] + xi[d] - xi[a] - xi[b];
            if *xi != base_xi {
                let new_xi: Vec<i64> = xi.iter().enumerate().map(|(i,&x)| diff_b * x - diff * base_xi[i]).collect();
                new_consv.push(new_xi);
            }
        }

        if self.is_separable(new_consv.clone()) {
            return false;
        }

        self.consv = new_consv.clone();
        return true;

    }

    fn get_edges(&mut self) -> Vec<Vec<(usize,usize)>> {
        let mut hm: HashMap<Vec<i64>,Vec<(usize,usize)>> = HashMap::new();

        for v in (0..self.n).combinations_with_replacement(2) {
            let mut consv_vector: Vec<i64> = vec![];
            for xi in &self.consv {
                consv_vector.push(xi[v[0]] + xi[v[1]]);
            }

            match hm.get_mut(&consv_vector) {
                Some(x) => {
                    x.push((v[0],v[1]));
                },
                None => {
                    hm.insert(consv_vector.clone(),vec![(v[0],v[1])]);
                }
            }
        }

        let mut edges:Vec<Vec<(usize,usize)>> = vec![];

        for (_, x) in hm.iter() {
            if x.len() > 1 {
                let mut x_sort = x.clone();
                x_sort.sort();
                edges.push(x_sort.clone());
            }
        }

        edges.sort();
        self.edges = edges.clone();

        return edges;
    }

}


struct InteractionsModEquiv {
    n: usize,
    my_inter_hs: HashSet<Vec<Vec<(usize,usize)>>>,
    my_inter_list: Vec<Interaction>,
    new_edge_list: Vec<(usize,usize,usize,usize)>,
}

impl InteractionsModEquiv {

    fn new(n:usize) -> Self {

        let mut new_edge_list: Vec<(usize,usize,usize,usize)> = vec![];
        for origin in (0..n).combinations_with_replacement(2) {
            let min = origin[0];
            for target in (min+1..n).combinations_with_replacement(2) {
                if origin[1] != target[0] && origin[1] != target[1] {
                    new_edge_list.push((origin[0],origin[1],target[0],target[1]));
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

    fn output_json(&mut self, file_name: String) -> std::io::Result<()> {
        
        // let serialized: String = serde_json::to_string_pretty(&(self.my_inter_list)).unwrap();
        let serialized: String = serde_json::to_string(&(self.my_inter_list)).unwrap();
    
        let mut file = File::create(file_name)?;
        file.write_all(serialized.as_bytes())?;
        Ok(())
    }

    fn edges_list(&mut self, consv: Vec<Vec<i64>>) -> Vec<Vec<(usize,usize)>> {
        let mut hm: HashMap<Vec<i64>,Vec<(usize,usize)>> = HashMap::new();
        for v in (0..self.n).combinations_with_replacement(2) {
            let mut temp_vector: Vec<i64> = vec![];
            for xi in consv.iter() {
                temp_vector.push(xi[v[0]] + xi[v[1]]);
            }

            match hm.get_mut(&temp_vector) {
                Some(x) => {
                    x.push((v[0],v[1]));
                },
                None => {
                    hm.insert(temp_vector.clone(),vec![(v[0],v[1])]);
                }
            }
        }

        let mut edges:Vec<Vec<(usize,usize)>> = vec![];

        for (_, x) in hm.iter() {
            if x.len() > 1 {
                let mut x_sort = x.clone();
                x_sort.sort();
                edges.push(x_sort.clone());
            }
        }

        edges.sort();
        return edges;

    }

    fn add(&mut self, mut new_inter: Interaction) -> bool {

        let edges = new_inter.get_edges();
        let consv = new_inter.consv.clone();
        
        for perm in (0..self.n).permutations(self.n) {
            let mut perm_consv: Vec<Vec<i64>> = vec![];
            for xi in consv.iter() {
                let mut perm_xi:Vec<i64> = vec![];
                for &x in perm.iter() {
                    perm_xi.push(xi[x]);
                }
                perm_consv.push(perm_xi.clone());
            }
            
            let edges_permuted = self.edges_list(perm_consv);

            match self.my_inter_hs.contains(&edges_permuted) {
                true => {
                    return false;
                },
                false => {}
            }
        }
        
        self.my_inter_list.push(new_inter);
        self.my_inter_hs.insert(edges.clone());
        return true;

    }

    fn create_list(&mut self) {
        let trivial_inter = Interaction::new(self.n);
        self.add(trivial_inter.clone());
        self.add_to_list(trivial_inter, 0);
    }

    fn add_to_list(&mut self, inter: Interaction, index: usize) {
        
        for i in index..self.new_edge_list.len() {
            let (a,b,c,d) = self.new_edge_list[i];
            let mut new_inter = inter.clone();
            if new_inter.merge((a,b),(c,d)) {
                if self.add(new_inter.clone()) {
                    self.add_to_list(new_inter.clone(), i+1);
                }
            }
        }
    }

}

fn main() {
    input!{
        n: usize,
    }

    //let mut inter = Interaction::new(n);
    let mut inter_list = InteractionsModEquiv::new(n);
    inter_list.create_list();
    let file_name = String::from("output.json");
    let result = inter_list.output_json(file_name);
    match result {
        Ok(..) => {},
        Err(err) => {
            println!("{err}");
        }
    }

}
