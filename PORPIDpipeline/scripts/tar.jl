using Pkg
Pkg.activate("/fh/fast/cohn_l/grp/Pipelines/PORPIDpipeline")
Pkg.instantiate()
using PORPIDpipeline, NextGenSeqUtils, FASTX

# zip porpid and postproc directories for easy download

function my_read_fasta_records(filename)
    stream = open(FASTA.Reader, filename)
    records = FASTA.Record[]
    for entry in stream
        push!(records, entry)
    end
    return records
end

function read_fasta_with_names_and_descriptions(filename; seqtype=String)
    records = my_read_fasta_records(filename)
    return [FASTA.identifier(r) for r in records],
            [FASTA.description(r) for r in records],
             seqtype[FASTA.sequence(seqtype, r) for r in records]
end

function write_fasta_with_names_and_descripts(filename::String, seqs; names = String[], descripts = String[])
    if length(names) > 0 && (length(names) != length(seqs) || length(names) != length(descripts) )
        error("number of sequences does not match number of names")
    end
    if length(names) == 0
        names = ["seq_$i" for i in 1:length(seqs)]
        descripts = ["" for i in 1:length(seqs)]
    end
    stream = open(FASTA.Writer, filename)
    for (name, seq, descript) in zip(names, seqs, descripts)
        write(stream, FASTA.Record(name, descript, seq))
    end
    close(stream)
end

dataset = snakemake.wildcards["dataset"]
porpid_dir = "porpid/$(dataset)"
postproc_dir = "postproc/$(dataset)"

degap_flag = snakemake.params["degap"]
println( "degap flag = $(degap_flag)")
samples = snakemake.params["samples"]
if degap_flag == "true"
    println("generating degapped postproc sequences...")
    run(`mkdir -p $(postproc_dir)/degapped_fasta/`)
    for sample in samples
        infile = "$(postproc_dir)/$(sample)/$(sample).fasta"
        outfile = "$(postproc_dir)/degapped_fasta/$(sample)-degapped.fasta"
        print(".")
        names, seqs = read_fasta_with_names(infile)
        names, descripts, seqs = read_fasta_with_names_and_descriptions(infile)
        write_fasta_with_names_and_descripts(outfile, degap.(seqs), names = names, descripts = descripts )
    end
    println()
end

# and now tar and zip both porpid and postproc directories
# first copy some reports from porpid to postproc
println("archiving porpid ...")
run(`cp $(porpid_dir)/contam_report.csv $(postproc_dir)/contam_report.csv`)

# now rename porpid_dir, zip and rename back
run(`mv $(porpid_dir) $(porpid_dir)-porpid`)
run(`tar -C porpid -czf porpid/$(dataset)-porpid.tar.gz $(dataset)-porpid`)
run(`mv $(porpid_dir)-porpid $(porpid_dir)`)
println("porpid directory archived and zipped ...")

# now do the same with postproc_dir
println("archiving postproc ...")
run(`mv $(postproc_dir) $(postproc_dir)-postproc`)
run(`tar -C postproc -czf postproc/$(dataset)-postproc.tar.gz $(dataset)-postproc`)
run(`mv $(postproc_dir)-postproc $(postproc_dir)`)
println("postproc directory archived and zipped ...")
