process mdust {
    tag "mdust-${file}"
    input: 
        path(file)
    output:
        path("${file}.dusted")
    conda 'bioconda::mdust'

"""
mdust ${file} > ${file}.dusted
"""
}

// TODO: Everytime this is referenced, see if any of the parameters need to be a variable
// and update accordingly.
process cleanup_tandem {
    tag "cleanup_tandem-${file}"
    input:
        path(file)
    output:
        path("${file}.cleaned")
    script:
"""
perl ${projectDir}/util/cleanup_tandem.pl \
    -misschar N \
    -nc 50000 \
    -nr 0.9 \
    -minlen 100 \
    -minscore 3000 \
    -trf 1 \
    -cleanN 1 \
    -cleanT 1 \
    -f ${file} > ${file}.cleaned

${file} > ${file}.pass.fa
"""
}

process tesorter {
    tag "tesorter-${file}"
    input:
        path(file)
    output:
        path("${file}.tesorted")
    conda 'bioconda::tesorter'
    cpus 4
"""
TEsorter ${file} --disable-pass2 -p ${task.cpus} > ${file}.tesorted
"""
}

process cleanup_misclassified {
    tag "cleanup_misclassified-${file}"
    input:
        path(file)
    output:
        path("${file}.cleaned")
    script:
// TODO: What is the output? What are the proper inputs?
"""
perl ${projectDir}/util/cleanup_misclas.pl ${file}
"""
}

process bed2gff {
    tag "bed2gff-${file}"
    input:
        path(file)
        val(type)
    output:
        path("${file.baseName}.gff3")
    script:
"""
perl ${projectDir}/util/bed2gff.pl ${file} ${type} > ${file.baseName}.gff3
"""
}