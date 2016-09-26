with open("testdata.txt", 'w') as outputfile:
    for x in range(2000):
        outputfile.write(str(x) + '\n')
