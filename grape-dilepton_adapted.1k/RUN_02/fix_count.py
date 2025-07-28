input_file = 'Grape.hepmc'
output_file = 'Grape_fixed_vertices.hepmc'

with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
    for line in fin:
        if line.startswith('E'):
            parts = line.strip().split()
            try:
                parts[2] = str(int(parts[2]) - 1)  # Decrement the middle number
            except (IndexError, ValueError):
                print("Warning: could not parse line:", line.strip())
            fout.write(' '.join(parts) + '\n')
        else:
            fout.write(line)
