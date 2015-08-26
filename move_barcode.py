#Function to move 8nt barcode from reverse reads to forward reads

def move_barcode(infile_1, infile_2, directory):
    with open(infile_1) as in1:
        with open (infile_2) as in2:
            with open(directory + 'R1.fastq', 'w') as out1:
                with open(directory + 'R2.fastq','w') as out2:
                    for line1 in in1:
                        line2 = in2.readline()
                        if line1[1:4] == 'SRR':
                            out1.write(line1)
                            out2.write(line2)
                        else:
                            out1.write(line2[:8] + line1)
                            out2.write(line2[8:])


