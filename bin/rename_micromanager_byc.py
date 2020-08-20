from byc.file_management import rename_byc
import sys

def run():

    if len(sys.argv) > 1:
        print(sys.argv)
        assert len(sys.argv) == 2, "Enter one directory as arg variable"
        expt_dir = sys.argv[1]
        rename_byc(expt_dir=expt_dir)
    else:
        # Don't pass rename_byc() an expt_dir keyword arg,
        # instead rename_byc() will have them choose it manually
        # in tkinter dialog
        rename_byc()

if __name__ == "__main__":
    run()