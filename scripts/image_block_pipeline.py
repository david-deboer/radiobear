import planet
j = planet.Planet('Jupiter')
frq = 43.34
bstep = 0.01
number_of_blocks = 15
chunking = 4
block_no = int(raw_input("Block number (N of {:.0f}): ".format(int(number_of_blocks / 4) + 1)))
thisrun = range((block_no - 1) * chunking + 1, block_no * chunking + 1)

for i in thisrun:
    if i > number_of_blocks:
        break
    j.run(freqs=frq, b=bstep, block=[i, number_of_blocks])
