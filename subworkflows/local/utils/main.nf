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