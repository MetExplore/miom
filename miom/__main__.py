import argparse
import miom

def convert_gem(args):
    print(f"MIOM v{miom.__version__}")
    m = miom.mio.load_gem(args.input)
    print(f"Loaded network with {m.num_reactions} reactions")
    print(f"Exporting to {args.output}...")
    miom.mio.export_gem(m, args.output)
    print("Done.")

def get_args():
    parser = argparse.ArgumentParser(description=f"MIOM: Mixed Integer Optimization for Metabolism")
    parser.add_argument(
        "--version",
        action="version", 
        version=f"MIOM v{miom.__version__}"
    )
    subparsers = parser.add_subparsers(help="sub-command help")
    convert = subparsers.add_parser("convert", help="Convert a model to miom format")
    convert.add_argument(
        "input",
        help="Input model file (if cobra is installed, any format supported by cobra is allowed)"
    )
    convert.add_argument(
        "output",
        help="Output GEM file in MIOM format"
    )
    convert.set_defaults(func=convert_gem)
    return parser.parse_args()

if __name__ == '__main__':
    args = get_args()
    if args.func:
        args.func(args)