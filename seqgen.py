from libprep.dnagenerator import generate_dna
import numpy as nu
import yaml

# load parameters
with open(r'generator.yaml') as file:
    parameters = yaml.load(file, Loader=yaml.FullLoader)
print("\nParameters:\n")
print(yaml.dump(parameters))
dna = generate_dna(parameters)

print("Sequence genrated!")