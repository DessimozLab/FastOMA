import argparse
from ._wrappers import logger
from .zoo.utils import auto_open


def extract_pw_rels(args):
    from lxml import etree
    from .zoo.hog import transform
    xml = etree.parse(args.orthoxml)
    with auto_open(args.out, 'wt') as fout:
        for p1, p2 in transform.iter_pairwise_relations(xml, rel_type=args.type, id_attribute="protId"):
            fout.write(f"{p1}\t{p2}\n")


def main():
    parser = argparse.ArgumentParser(description="FastOMA helper scripts")
    parser.add_argument('-v', default=0, action="count", help="increase verbosity")
    subparsers = parser.add_subparsers(required=True)

    parser_pw = subparsers.add_parser('pw-rel')
    parser_pw.add_argument("--type", choices=("ortholog", "paralog"), default="ortholog",
                           help="Type of relations to extract. either 'ortholog' or 'paralog'")
    parser_pw.add_argument("--out", required=True, help="Path to output file")
    parser_pw.add_argument("--orthoxml", required=True, help="Path to input orthoxml file")
    parser_pw.set_defaults(func=extract_pw_rels)

    conf = parser.parse_args()
    logger.setLevel(level=30 - 10 * min(conf.v, 2))
    logger.debug(conf)
    conf.func(conf)


if __name__ == "__main__":
    main()