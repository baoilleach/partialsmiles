from smiparser import ParseSmiles

def partialtest():
    for line in smis.split("\n"):
        print(line)
        for i in range(1, len(line)):
            try:
                mol = ParseSmiles(line[:i], partial=True)
            except Exception as e:
                print(str(e))
                break

        try:
            mol = ParseSmiles(line, partial=False)
        except Exception:
            pass
        else:
            print("There should have been an exception")

def main():
    with open(r"C:\Tools\LargeData\sortedbylength.smi") as inp:
        for line in inp:
            smi = line.split()[0]
            try:
                mol = ParseSmiles(smi, partial=False)
            except Exception as e:
                print(str(e))

if __name__ == "__main__":
    # partialtest()
    main()
