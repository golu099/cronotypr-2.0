version 1.0 

import "../tasks/tasks_tools.wdl" as tools

workflow plasmidotyper {
    meta {
        description: " Running plasmidotyper change later description."
    }
    call tools.blastn_db {
        input:
            database_fastas = database_fastas,
            query_fastas = query_fastas
    }
    call tools.blast_to_excel {
        input: 
            blast_file = blast_file
    }
    call tools.genotyping_length {
        input: 
            blast_csv = blast_to_excel.blast_csv


    }
}
