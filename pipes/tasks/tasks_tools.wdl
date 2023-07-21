version 1.0 
#make db of plasmids and run  blastn 
task blastn_db {
    input {
        Array[File] database_fastas
        Array[File] query_fastas
        Int     perc_identity = 50
        Int     best_hit_overhang = 0.25
        Int     best_hit_score_edge = 0.05
        String     evalue = 1e-10
        Int     num_threads = 12
        Int?    machine_mem_gb
        Int     disk_size_gb = ceil(2*size(otu_ref, "GiB")) + 5
        String  docker = "quay.io/broadinstitute/ncbi-tools:2.10.7.10"
        #what docker image should i use? Docker should contain ncbi and python 
        #String  docker 
    }
    meta {
        description: "Make a database from multiple files, and then run them through blastn"
    }
    command <<<
        set -ex -o pipefail
        #make new list of database fasta names 
        query_names_list= list_sample_names.txt
        for x in ~{sep=' ' query_fastas}; do
            NEWSAMPLENAME=$(basename $fasta .fasta | perl -lape 's/[_]/-/g')
            echo $NEWSAMPLENAME >> list_sample_names.txt
        #concat 
        cat ~{sep=" " database_fastas} > concat_db.fasta
        cat ~{sep=" " query_fastas} > concat_query.fasta
        #add the concat file as the "nt.fsa"
        blastn –db nt –query concat_db.fasta –out concat_db_ncbi
        blastn -query concat_query.fasta -out blast_output.txt -task megablast -db concat_db_ncbi ~{"-num_threads" + num_threads} ~{"-evalue" + evalue} ~{"-best_hit_score_edge" + best_hit_score_edge} ~{"-best_hit_overhang"+ best_hit_overhang} -outfmt 6 ~{"-perc_identity" + perc_identity}
>>>

output { 
    File    raw_blast = blast_output.txt
    File    query_names_list = list_sample_names.txt
}
runtime {
    docker: docker
    cpu: 2
    memory: select_first([machine_mem_gb, 6]) + "GB"
    disk: disk_size_gb + " GB"
    disks: "local-disk " + disk_size_gb + " HDD"
}

}
task blast_to_excel{
    input {
        File    raw_blast
        Int     cpu = 5
        Int     memory_mb = 2000 
        Int     disk_size_gb = ceil(2*size(raw_blast, "GiB")) + 5  
        String  
    }
    meta {
        description: "Takes in blast_results.tsv and converts it into a .csv file."
    }
    parameter_meta {
        raw_blast: {
            description: "Import your blast file in .tsv format",
            category: "required"
        }
    }
    command <<<
        python3 << CODE 
        import sys
        import os 
        import re
        import csv
        with open(~{raw_blast}, "r") as text_file: 
            csv_reader: csv.reader(text_file, delimiter="\t")
            lines = text_file.readlines()[0:]
            mylist = []
            mylist.append("Strain	ST	Species_Strain	Species	Type	Plasmid	Gene	%ID	Hit_Length	Mismatches	Gaps	q.start	q.end	s.start	s.end	e-value	bit")
            for i in lines: 
                if '.' in i: 
                mylist.append(i)
        print("Performing conversion for "+ str(~{raw_blast})+ "\n")
        wr = csv.writer("$(date --iso-8601)_plasmidotyoer.csv", delimiter = ",")
        wr.writerows([elem.replace('\t','|') for elem in mylist])
    text_file.close()
    CODE    
>>>
    output {
        File    blast_csv = "$(date --iso-8601)_blast_output.csv"
    }
    runtime { 
        docker: docker
        memory: "${memory_mb} MiB"
        cpu: cpu
        disk: disk_size_gb + "GB"
        disks: "local-disk" + disk_size_gb + "HDD"
    } 
}
task genotyping_length {
    input{
        File    blast_csv
        #in the future have a list of names for each sample made already in the first step so that it doesnt have to be uploaded as a separate name file 
        File    query_names_list
        String  
    }
    meta {
        description: "Moving data into a plasmidotyper table."
    }
    command <<<
        python3<<< CODE
        import sys
        import os 
        import re
        import csv
        import pandas as pd
        import numpy as np
        # argument 1= excel file
        # argument 2= modified excel output
        # argument 3=list.txt
        mylist = []
        with open(~{blast_csv}, "r") as text_file, open(~{query_names_list}, "w+") as outfile_2:
            data = pd.read_csv(~{blast_csv}, header=0, sep=",",
                                error_bad_lines=False, engine='python')
        #Change the header of the query_fastas to inlcude only column name. 
        #Then made a list of "Stain Names" to dataframe and deleted unecessary columns
            for x in data['Strain']:
            b = (re.search(r'(\w.+\d)+(_+\w*)', str(x)))
            if b:
                outfile_2.write(str(b.group(1))+'\n')
            else:
                p = (re.search(r'(\w.+\d.)+(_+\w*)', str(x)))
                if p:
                    outfile_2.write(str(p.group(1))+'\n')

                else:
                    l = (re.search(r'(\w.+\d.)+(c+\w*)', str(x)))
                    if l:
                        outfile_2.write(str(l.group(1))+'\n')
                    else:
                        outfile_2.write(str(x)+'\n')
        # remove blanke lines from list
    with open("$(date --iso-8601)_plasmidotyper.csv", "wb") as bank, open(~{query_names_list}, "r+") as list1:
        list1 = list1.readlines()
        data['Strain_name'] = pd.Series(list1)
    #finish adding list to dataframe 
        del data['Strain']
        del data['Mismatches']
        del data['Gaps']
        del data['s.start']
        del data['s.end']
        del data['e-value']
        del data['bit']
        print('Computing results for ' + "~{blast_csv}")
    #data_pivot_table 
    #set max ID to 100% 
        data = data.pivot_table(index=['Species', 'Type', 'Plasmid', 'Gene', 'Hit_Length', 'q.start', 'q.end'], columns=[
                            'Species_Strain', 'ST', 'Strain_name'], values=['%ID'], aggfunc='sum', fill_value=0)
    data = data.clip(upper=100)
    data.to_csv(sys.argv[2], index=True)
    outfile_2.close()
    text_file.close()
    outfile.close()
    CODE
>>>
    output {
        File    plasmidotyper_results = "$(date --iso-8601)_plasmidotyper.csv"    
}
    runtime {
        docker: docker
        memory: "~{memory_mb} GB"
        cpu: cpu
        disk: disk_size_gb + " GB"
        disks: "local-disk " + disk_size_gb + " HDD"
    }
}