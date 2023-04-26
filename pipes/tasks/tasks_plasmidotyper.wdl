version 1.0 
#argument 1= Blast txt file 
#argument 2 = Excel output file
task blast_to_excel{
    input {
        File    blast_file
        Int     cpu = 5
        Int     memory_mb = 2000 
        Int     disk_size_gb = ceil(2*size(blast_file, "GiB")) + 5  
    }
    parameter_meta {
        blast_file: {
            description: "Import your blast file in .tsv format",
            category: "required"
        }
    }
    command <<<
        python3 << CODE 
        with open(~{blast_file}, "r") as text_file: 
            csv_reader: csv.reader(text_file, delimiter="\t")
            lines = text_file.readlines()[0:]
            mylist = []
            mylist.append("Strain	ST	Species_Strain	Species	Type	Plasmid	Gene	%ID	Hit_Length	Mismatches	Gaps	q.start	q.end	s.start	s.end	e-value	bit")
            for i in lines: 
                if '.' in i: 
                mylist.append(i)
        print("Performing conversion for "+ str(~{blast_file})+ "\n")
        wr = csv.writer("$(date --iso-8601)_plasmidotyoer.csv", delimiter = ",")
        wr.writerows([elem.replace('\t','|') for elem in mylist])
    text_file.close()
    CODE    
>>>
    output {
        File    csv_file = "$(date --iso-8601)_plasmidotyoer.csv"
    }
    runtime { 
        memory: "${memory_mb} MiB"
        cpu: cpu
        disk: disk_size_gb + "GB"
        disks: "local-disk" + disk_size_gb + "HDD"
    } 
}