use proconio::input;
mod interaction;

fn main() {
    input! {
        n: usize,
    }

    let mut inter_list = interaction::InteractionsModEquiv::new(n);
    inter_list.create_list();
    let mut file_name = String::from("output/size_");
    file_name.push_str(&(n as i32).to_string());
    file_name.push_str(".json");
    let result = inter_list.output_json(file_name);
    match result {
        Ok(..) => {}
        Err(err) => {
            println!("{err}");
        }
    }
}
