#!/usr/bin/python3
import random,math
boxes = [['gold','gold'], ['gold','silver'],['silver','silver']]

N_TRIALS = 1000000
success = 0
fail = 0
total = 0

for trial in range(N_TRIALS):
    random_box = math.floor(random.random()*len(boxes))
    random_ball = math.floor(random.random()*len(boxes[random_box]))
    if boxes[random_box][random_ball] == 'silver':
        continue
    else:
        total = total + 1
        other_ball = (random_ball^1)
        if boxes[random_box][other_ball] == 'gold':
            success = success + 1
        else:
            fail = fail + 1

print(N_TRIALS," trials made")
print(total," times gold drawn first")
print(success," gold drawn second - ",round(success/total,4),"%")
print(fail," silver drawn second - ",round(fail/total,4),"%")