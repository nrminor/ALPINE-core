pub mod lib {

    pub fn find_input_path(args: &[String]) -> &String {

        let fasta_path = &args[1];
        fasta_path
    
    }

    pub mod distmat;
    pub mod filtering;
    pub mod prevalence;
    pub mod reporting;

}
