
// lib/Utils.groovy
import java.nio.file.Path

class Utils {
    static def getMaxFileSize(Path folderPath) {
        def maxFileSize = 0L
        def traverseFolder = { Path currentPath ->
            def currentFile = currentPath.toFile()
            if (currentFile.isDirectory()) {
                def files = currentFile.listFiles()
                if (files) {
                    files.each { file ->
                        if (file.isDirectory()) {
                            traverseFolder(file.toPath())
                        } else {
                            def fileSize = file.length()
                            if (fileSize > maxFileSize) {
                                maxFileSize = fileSize
                            }
                        }
                    }
                }
            } else {
                def fileSize = currentFile.length()
                if (fileSize > maxFileSize) {
                    maxFileSize = fileSize
                }
            }
        }
        traverseFolder(folderPath)
        return maxFileSize
    }

    static def mem_cat(filesize, nr_genomes) {
        def fac = 1
        if (nr_genomes > 500) fac = 2
        if (filesize < 1000000) return 12.GB * fac
        else if (filesize < 2000000) return 20.GB * fac
        else if (filesize < 5000000) return 32.GB * fac
        else return 64.GB * fac
    }

    static def time_cat(filesize, nr_genomes) {
        def fac = 1
        if (nr_genomes > 500) fac = 2
        if (filesize < 1000000) return 4.h * fac
        else if (filesize < 20000000) return 24.h * fac
        else return 72.h * fac
    }
}