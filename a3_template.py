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
                prior_prob = math.log(10,self.prior[states[i]])
                prob_emission = math.log(10,self.emission[states[i]][sequence[i]])
                prob_of_state = prior_prob * prob_emission
                probability.append(prob_of_state)
        # Start your code
            else:
                prob_emission = math.log(10,self.emission[states[i]][sequence[i]])
                #transit from previous state to current
                prob_transition = math.log(10,self.transition[states[i-1]][states[i]])
                prob_of_state = probability[i-1]*prob_transition*prob_emission
                probability.append(prob_of_state)
        log_prob =1
        for prob in probability:
            print(prob)
            log_prob = log_prob * prob
        print("My code here")
        return log_prob
        # End your code
        ###########################################


    # Outputs the most likely sequence of states given an emission sequence
    # - sequence: String with characters [A,C,T,G]
    # return: list of state indices, e.g. [0,0,0,1,1,0,0,...]
    def viterbi(self, sequence):
        ###########################################
        # Start your code
        output_list=[]
        for x in range(len(sequence)):
            #print(sequence[x])
            argmax = []
            for state in range(self.num_states):
                new_state =[]
                next_state = copy.deepcopy(output_list)
                next_state.append(state)
                argmax.append(self.logprob(sequence,next_state))
            #print(argmax.index(max(argmax)))
            output_list.append(argmax.index(max(argmax)))
           # print(output_list)
        print("My code here")
        return output_list
        # End your code
        ###########################################

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
logprob = hmm.logprob(sequence, viterbi)
print(logprob)
write_output("my_small_output.txt", logprob, viterbi)


#sequence = read_sequence("ecoli.txt")
#viterbi = hmm.viterbi(sequence)
#logprob = hmm.logprob(sequence, viterbi)
#write_output("ecoli_output.txt", logprob, viterbi)


