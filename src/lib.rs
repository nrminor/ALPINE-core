pub mod lib {

    pub fn find_input_path(args: &[String]) -> &String {
        &args[1] as _
    }

    pub mod distmat;
    pub mod filtering;
    pub mod prevalence;
    pub mod reporting;
}
