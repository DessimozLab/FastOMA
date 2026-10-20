// Order- and id-independent views of FastOMA output files for nf-test snapshots.
//
// HOG/OG numbers and the row order of the TSV files depend on task completion and set
// iteration order, so the snapshots compare group memberships, ortholog pairs and the
// nesting of the orthoXML groups instead of raw file content.

import java.nio.file.Path
import java.util.zip.GZIPInputStream
import javax.xml.parsers.DocumentBuilderFactory
import org.w3c.dom.Element
import org.w3c.dom.Node

class FastomaOutput {

    static List<String> readLines(Path file) {
        if (file.toString().endsWith('.gz')) {
            return new GZIPInputStream(file.newInputStream()).withReader { it.readLines() }
        }
        return file.readLines()
    }

    // Groups of a "group<TAB>member" TSV (RootHOGs.tsv, OrthologousGroups.tsv) as a sorted
    // list of sorted member lists; the group names are dropped.
    static List<List<String>> groupPartition(Path tsv) {
        def groups = [:].withDefault { [] as SortedSet }
        readLines(tsv).drop(1).findAll { it }.each { line ->
            def cols = line.split('\t')
            groups[cols[0]] << cols[1]
        }
        return groups.values().collect { it as List }.sort { it.join(',') }
    }

    // Unordered pairs of a two-column TSV (orthologs.tsv.gz) as sorted "a <-> b" strings.
    static List<String> pairs(Path tsv) {
        return readLines(tsv).findAll { it }.collect { line ->
            def p = line.split('\t')[0..1].sort()
            "${p[0]} <-> ${p[1]}".toString()
        }.sort()
    }

    private static Element parseOrthoxml(Path orthoxml) {
        def factory = DocumentBuilderFactory.newInstance()
        factory.setNamespaceAware(true)
        return factory.newDocumentBuilder().parse(orthoxml.toFile()).documentElement
    }

    private static List<Element> childElements(Node node) {
        def out = []
        def children = node.childNodes
        for (int i = 0; i < children.length; i++) {
            if (children.item(i) instanceof Element) {
                out << (Element) children.item(i)
            }
        }
        return out
    }

    private static Map<String, String> geneIdToProtId(Element root) {
        def map = [:]
        def genes = root.getElementsByTagNameNS('*', 'gene')
        for (int i = 0; i < genes.length; i++) {
            def g = (Element) genes.item(i)
            map[g.getAttribute('id')] = g.getAttribute('protId')
        }
        return map
    }

    private static String canonicalGroup(Element elem, Map<String, String> id2prot) {
        switch (elem.localName) {
            case 'geneRef':
                return id2prot.get(elem.getAttribute('id'), elem.getAttribute('id'))
            case 'orthologGroup':
            case 'paralogGroup':
                def taxRange = childElements(elem).find {
                    it.localName == 'property' && it.getAttribute('name') == 'TaxRange'
                }?.getAttribute('value')
                def children = childElements(elem)
                    .findAll { it.localName in ['orthologGroup', 'paralogGroup', 'geneRef'] }
                    .collect { canonicalGroup(it, id2prot) }
                    .sort()
                def label = elem.localName == 'orthologGroup' ? "OG[${taxRange}]" : (taxRange ? "PG[${taxRange}]" : 'PG')
                return "${label}(${children.join(', ')})".toString()
            default:
                return ''
        }
    }

    // One string per rootHOG describing the nesting of ortholog/paralog groups and their
    // genes (as protIds); group ids, scores and properties other than TaxRange are dropped.
    static List<String> canonicalHogs(Path orthoxml) {
        def root = parseOrthoxml(orthoxml)
        def id2prot = geneIdToProtId(root)
        def groups = childElements(root).find { it.localName == 'groups' }
        if (groups == null) {
            return []
        }
        return childElements(groups).collect { canonicalGroup(it, id2prot) }.sort()
    }

    // How often each protein is referenced from the groups; a protein must appear at most once.
    static Map<String, Integer> geneRefCounts(Path orthoxml) {
        def root = parseOrthoxml(orthoxml)
        def id2prot = geneIdToProtId(root)
        def counts = [:].withDefault { 0 }
        def refs = root.getElementsByTagNameNS('*', 'geneRef')
        for (int i = 0; i < refs.length; i++) {
            def id = ((Element) refs.item(i)).getAttribute('id')
            counts[id2prot.get(id, id)] += 1
        }
        return counts
    }

    static Map<String, Integer> elementCounts(Path orthoxml) {
        def root = parseOrthoxml(orthoxml)
        return ['orthologGroup', 'paralogGroup', 'gene', 'geneRef'].collectEntries { name ->
            [name, root.getElementsByTagNameNS('*', name).length]
        }
    }
}
