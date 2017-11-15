import random
import math
import copy

#####################################################
#####################################################
# Please enter the number of hours you spent on this
# assignment here
num_hours_i_spent_on_this_assignment = 0
#####################################################
#####################################################

#####################################################
#####################################################
# Give one short piece of feedback about the course so far. What
# have you found most interesting? Is there a topic that you had trouble
# understanding? Are there any changes that could improve the value of the
# course to you? (We will anonymize these before reading them.)
# <Your feedback goes here>
#####################################################
#####################################################



# Outputs a random integer, according to a multinomial
# distribution specified by probs.
def rand_multinomial(probs):
    # Make sure probs sum to 1
    assert(abs(sum(probs) - 1.0) < 1e-5)
    rand = random.random()
    for index, prob in enumerate(probs):
        if rand < prob:
            return index
        else:
            rand -= prob
    return 0

# Outputs a random key, according to a (key,prob)
# iterator. For a probability dictionary
# d = {"A": 0.9, "C": 0.1}
# call using rand_multinomial_iter(d.items())
def rand_multinomial_iter(iterator):
    rand = random.random()
    for key, prob in iterator:
        if rand < prob:
            return key
        else:
            rand -= prob
    return 0

class HMM():
    #0.5*0.169*0.01*0.169
    def __init__(self):
        self.num_states = 2
        self.prior = [0.5, 0.5]
        self.transition = [[0.999, 0.001], [0.01, 0.99]]
        self.emission = [{"A": 0.291, "T": 0.291, "C": 0.209, "G": 0.209},
                         {"A": 0.169, "T": 0.169, "C": 0.331, "G": 0.331}]
        self.path=[]
    # Generates a sequence of states and characters from
    # the HMM model.
    # - length: Length of output sequence
    def sample(self, length):
        sequence = []
        states = []
        rand = random.random()
        cur_state = rand_multinomial(self.prior)
        for i in range(length):
            states.append(cur_state)
            char = rand_multinomial_iter(self.emission[cur_state].items())
            sequence.append(char)
            cur_state = rand_multinomial(self.transition[cur_state])
        return sequence, states

    # Generates a emission sequence given a sequence of states
    def generate_sequence(self, states):
        sequence = []
        for state in states:
            char = rand_multinomial_iter(self.emission[state].items())
            sequence.append(char)
        return sequence

    # Computes the (natural) log probability of sequence given a sequence of states.
    def logprob(self, sequence, states):
        ###########################################
        probability =[]

        for i in range(len(states)):
            if i ==0:
                prior_prob = math.log(self.prior[states[i]])
                prob_emission = math.log(self.emission[states[i]][sequence[i]])
                prob_of_state = prior_prob * prob_emission
                probability.append(prob_of_state)
        # Start your code
            else:
                prob_emission = math.log(self.emission[states[i]][sequence[i]])
                #transit from previous state to current
                prob_transition = math.log(self.transition[states[i-1]][states[i]])
                prob_of_state = probability[i-1]*prob_transition*prob_emission
                probability.append(prob_of_state)
        log_prob =0
        for prob in probability:
            log_prob = log_prob + prob
        print("My code here")
        return log_prob
        # End your code
        ###########################################


    # Outputs the most likely sequence of states given an emission sequence
    # - sequence: String with characters [A,C,T,G]
    # return: list of state indices, e.g. [0,0,0,1,1,0,0,...]
    def viterbi(self, sequence):
        ###########################################
        # Start your cod
        steps = [][]
        steps[0][0]= self.prior[0] *self.emission[0][sequence[0]]
        steps[0][1] = self.prior[1] *self.emission[1][sequence[0]]
        for sym in range(1,len(sequence)):
            # 1= move, 0 = stay
            L_to_L = steps[sym-1][0] * self.transition[0][0] *self.emission[0][sequence[sym]]
            H_to_L = steps[sym-1][1] * self.transition[1][1] *self.emission[0][sequence[sym]]
            H_to_H = steps[sym-1][1] * self.transition[1][0] *self.emission[1][sequence[sym]]
            L_to_H = steps[sym-1][0] * self.transition[0][1] *self.emission[1][sequence[sym]]
            steps[sym][0] = L_to_L + H_to_L
            steps[sym][1] = L_to_H +H_to_H
        return steps

                # need to consider transition probabilities
                # probability of state 0 from previous probabilities

                #probability of state 0 give symb[i] when intial state came from H (H->L)

                #probability of state 0 givent sym[i] when intial state came from L (H->L)
                 # SUM
                #probability of state 1 given symb[i] when intial state came from H (H->L)

                #probability of state 1 given symb[i] when intiatle state came from L (L->L)
                    #SUM
                # store max

        # End your code
        ###########################################
    def getSequence(self,path):
        outputsequence=""
        for tup in path:
            outputsequence = outputsequence + tup[0]
            print(outputsequence)
        return outputsequence
    def obtain_max(self,probL, probR):
        if max(probL, probR) == probL:
            return ('0', probL)
        else:
            return ('1', probR)

def read_sequence(filename):
    with open(filename, "r") as f:
        return f.read().strip()

def write_sequence(filename, sequence):
    with open(filename, "w") as f:
        f.write("".join(sequence))

def write_output(filename, logprob, states):
    with open(filename, "w") as f:
        f.write(str(logprob))
        f.write("\n")
        for state in range(2):
            f.write(str(states.count(state)))
            f.write("\n")
        f.write("".join(map(str, states)))
        f.write("\n")

hmm = HMM()

sequence = read_sequence("small.txt")
viterbi = hmm.viterbi(sequence)
print(viterbi)
#logprob = hmm.logprob(sequence, viterbi)

write_output("my_small_output.txt", 0, viterbi)


#sequence = read_sequence("ecoli.txt")
#viterbi = hmm.viterbi(sequence)
#logprob = hmm.logprob(sequence, viterbi)
#write_output("ecoli_output.txt", logprob, viterbi)


