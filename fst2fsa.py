# -*- coding: utf-8 -*-
from __future__ import print_function
import sys
import argparse
from collections import defaultdict

def print_state(state, transitions, separator):
    print(state, "", 0, sep=separator)
    print(state, end="")
    for t in transitions:
        for i in range(len(transitions[t])):
            state_name = state+"_"+t+"_"+str(i+1)
            print(separator + state_name + separator+"0", end="")
    print("")
    for t in transitions:
        for i in range(len(transitions[t])):
            state_name = state+"_"+t+"_"+str(i+1)
            print(state_name, transitions[t][i], 0, sep=separator)
            print(state_name, t, 0, sep=separator)

def main(args):
    print(args.separator)
    print(args.start)
    print(args.end)
    
    print(args.start, "", 0, sep=args.separator)
    print(args.start, "0", 0, sep=args.separator)
    
    def is_epsilon(s):
        return len(s) > 2 and s[0] == s[-1] and s[0] == "@"
    
    transitions = defaultdict(lambda: defaultdict(list))
    
    for line in sys.stdin:
        line = line.strip().split(args.separator)
        # if line[0] != prev_state:
            # if prev_state != "":
                # print_state(prev_state, transitions, args.separator)
            # prev_state = line[0]
            # transitions = defaultdict(list)
        if len(line) < 4:
            # final state
            transitions[line[0]][args.end].append("")
        else:
            transitions[line[0]][line[1]].append("" if is_epsilon(line[2]) else line[2])
    for state in transitions:
        print_state(state, transitions[state], args.separator)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
                description="""Author: Gábor Borbély, License: MIT
Contact: http://math.bme.hu/~borbely/indexeng""",
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("-s", "--separator", dest='separator', type=str, default="\t",
                        help="")
    parser.add_argument("--start", dest='start', type=str, default="^", help="")
    parser.add_argument("--end", dest='end', type=str, default="$", help="")
    # parser.add_argument("--epsilon", dest='epsilon', type=str, default=["@0@"], help="", nargs="+")
                
    exit(main(parser.parse_args()))
