import argparse


def some_function(target, end="!"):
    """Some example funcion"""
    msg = "hi " + target + end
    print(msg)


def start():
    # All the logic of argparse goes in this function
    parser = argparse.ArgumentParser(description="Say hi.")
    parser.add_argument("target", type=str, help="the name of the target")
    parser.add_argument(
        "--end",
        dest="end",
        default="!",
        help="sum the integers (default: find the max)",
    )

    args = parser.parse_args()
    some_function(args.target, end=args.end)
