// Test-only shim: nf-test's `nextflow_function` test type generates an
// `include { <name> } from '<script>'`, but Nextflow's module loader only
// exposes top-level `def`/`process`/`workflow` declarations as includable
// components -- a bare `class Utils { static ... }` (as in ../../lib/Utils.groovy)
// cannot be included directly. These thin top-level wrappers make Utils'
// static methods includable for lib/Utils.groovy.test without changing the
// production code.

def mem_cat(filesize, nr_genomes) {
    return Utils.mem_cat(filesize, nr_genomes)
}

def time_cat(filesize, nr_genomes) {
    return Utils.time_cat(filesize, nr_genomes)
}

def getMaxFileSize(folderPath) {
    return Utils.getMaxFileSize(folderPath)
}
