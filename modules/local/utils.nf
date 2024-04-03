def idFromFileName(fileName) {

    def trial = ( fileName
        ).replaceFirst(
            /\.f(ast)?q$/, ''
        ).replaceFirst(
            /\.f(asta|sa|a|as|aa)?$/, ''
        ).replaceFirst(
            /\.gff(3)?$/, ''
        ).replaceFirst(
            /\.gz$/, ''
        )

    if ( trial == fileName ) { return fileName }

    return idFromFileName ( trial )
}