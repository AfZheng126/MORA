
/*
Sorts the Vector in terms of the first element in the Tuple
We are assuming tha the first elements in the tuples are all unique from each other
 */
pub(crate) fn sort_vec(vec: &Vec<(usize, f32)>) -> Vec<(usize, f32)> {
    let mut vec = vec.to_owned();
    let mut sorted_list = Vec::new();
    while !vec.is_empty() {
        let mut min = vec[0];
        for i in 0..vec.len() {
            if vec[i].0 < min.0 {
                min = vec[i];
            }
        }
        sorted_list.push(min);
        let index = vec.iter().position(|&x| x.0 == min.0).unwrap();
        vec.remove(index);
    }
    sorted_list
}
