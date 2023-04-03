import random
import sys

n = int(sys.argv[1])

with open("input.txt", "w") as f:
    f.write(str(n) + "\n")
    for i in range(n):
        random_int = random.randint(0, 9)
        if i == n-1:
            f.write(str(random_int))
        else:
            f.write(str(random_int) + "\n")

print(f"{n} random numbers have been generated and saved to input.txt")
