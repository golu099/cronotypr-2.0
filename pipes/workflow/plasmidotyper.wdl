version 1.0 

import "../tasks/tasks_tools.wdl" as tools

workflow plasmidotyper {
    meta {
        description: " Running plasmidotyper change later description."
    }

    call tools.blast_to_excel {
        input: 
            blast_file = blast_file
    }
}
