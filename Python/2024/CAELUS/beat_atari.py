# Imports
import gym
import numpy as np
import tensorflow as tf
import random
from numpy import random
import keras
from keras import backend as K
from tqdm import tqdm
from PIL import Image

# Global variables / constants
N_ACTIONS = 4
SEED = 0

# Class definitions

# Class to store memory for experience replay concept
class RingBuf:
    def __init__(self, size):
        self.data = [None] * (size + 1)
        self.start = 0
        self.end = 0
        
    def append(self, element):
        self.data[self.end] = element
        self.end = (self.end + 1) % len(self.data)
        # end == start and yet we just added one element. This means the buffer has one too many element. Remove the first element by incrementing start.
        if self.end == self.start:
            self.start = (self.start + 1) % len(self.data)
        
    def __getitem__(self, idx):
        return self.data[(self.start + idx) % len(self.data)]
    
    def __len__(self):
        if self.end < self.start:
            return self.end + len(self.data) - self.start
        else:
            return self.end - self.start
        
    def __iter__(self):
        for i in range(len(self)):
            yield self[i]

    def sample(self, num):
        out = []
        for i in range(num):
            out.append(self.__getitem__(np.random.randint(0, self.__len__())))
        return out

# Function Definitions

# Note: pass in_keras=False to use this function with raw numbers of numpy arrays for testing
def huber_loss(a, b, in_keras=True):
    error = a - b
    quadratic_term = error*error / 2
    linear_term = abs(error) - 1/2
    use_linear_term = (abs(error) > 1.0)
    if in_keras:
        # Keras won't let us multiply floats by booleans, so we explicitly cast the booleans to floats
        use_linear_term = K.cast(use_linear_term, 'float32')
    return use_linear_term * linear_term + (1-use_linear_term) * quadratic_term

# Training the DQN model, simulates all actions for state and choses best rewards, note the various actions are passed through a mask to find the best Q value for an action in one iteration. 
def fit_batch(model, gamma, start_states, actions, rewards, next_states, is_terminal):
    """Do one deep Q learning iteration.
    
    Params:
    - model: The DQN
    - gamma: Discount factor (should be 0.99)
    - start_states: numpy array of starting states
    - actions: numpy array of one-hot encoded actions corresponding to the start states
    - rewards: numpy array of rewards corresponding to the start states and actions
    - next_states: numpy array of the resulting states corresponding to the start states and actions
    - is_terminal: numpy boolean array of whether the resulting state is terminal
    
    """

    # First, predict the Q values of the next states. Note how we are passing ones as the mask.
    next_Q_values = model.predict([next_states, np.ones(actions.shape)],batch_size = 32)
    # The Q values of the terminal states is 0 by definition, so override them
    next_Q_values[is_terminal] = 0
    # The Q values of each start state is the reward + gamma * the max next state Q value
    Q_values = rewards + gamma * np.max(next_Q_values, axis=1)
    # Fit the keras model. Note how we are passing the actions as the mask and multiplying
    # the targets by the actions.
    model.fit([start_states, actions], actions * Q_values[:, None], epochs=1, batch_size=len(start_states), verbose=0)

# Creating the DQN model neural network 
def atari_model(n_actions):
    
    # We assume a theano backend here, so the "channels" are first, changed to TF backend 
    ATARI_SHAPE = (210, 160, 3)

    # With the functional API we need to define the inputs.
    frames_input = keras.layers.Input(ATARI_SHAPE, name='frames')
    actions_input = keras.layers.Input((n_actions,), name='mask')

    # Assuming that the input frames are still encoded from 0 to 255. Transforming to [0, 1].
    normalized = keras.layers.Lambda(lambda x: x / 255.0)(frames_input)
    
    # "The first hidden layer convolves 16 8×8 filters with stride 4 with the input image and applies a rectifier nonlinearity."
    conv_1 = keras.layers.Conv2D(filters=16, kernel_size=(8, 8), strides=(4, 4), activation='relu')(normalized)
    # "The second hidden layer convolves 32 4×4 filters with stride 2, again followed by a rectifier nonlinearity."
    conv_2 = keras.layers.Conv2D(filters=32, kernel_size=(4, 4), strides=(2, 2), activation='relu')(conv_1)
    # Flattening the second convolutional layer.
    conv_flattened = keras.layers.Flatten()(conv_2)
    # "The final hidden layer is fully-connected and consists of 256 rectifier units."
    hidden = keras.layers.Dense(256, activation='relu')(conv_flattened)
    # "The output layer is a fully-connected linear layer with a single output for each valid action."
    output = keras.layers.Dense(n_actions)(hidden)
    # Finally, we multiply the output by the mask! ---------Replacing old merge() function
    filtered_output = keras.layers.Multiply()([output, actions_input])

    model = keras.Model(inputs=[frames_input, actions_input], outputs = filtered_output)
    optimizer = optimizer = keras.optimizers.RMSprop(learning_rate=0.00025, rho=0.95, epsilon=0.01)
    model.compile(optimizer, loss=huber_loss) # Can cahnge to Huber Loss
    return model

def choose_best_action(model, state):
    # finds the maximum model predict value given state and list of predifined actions
    
    # create an instance of action space and state space
    action_space = []
    for i in range(0,N_ACTIONS):
        action_space.append(i)
    state = np.array([state])
    # call the model to predict Q-values for all actions
    Q_values = model.predict([state, action_space], batch_size = 1)
    # return array of Q-values and actions, conduct argmax for action
    maxIndex = np.where(Q_values == np.max(Q_values))
    return maxIndex[1][0]

def get_epsilon_for_iteration(iteration):
    # we want the rate at which the model chooses a random function to decrease over time to prioritize its more accurate learning
    start = 1 # starts at one
    decrease_rate = 9*(10**(-7))
    epsilon = start - decrease_rate*iteration
    if epsilon < .1: # epsilon will have a minimum value of .1
        epsilon = .1
    return epsilon

# Experience replay, at each step, play one step in the game and add all taken state from step just taken to table, then fit based on the memory, this table includes the rewards for each action, i.e. memory instance. 
def q_iteration(env, model, state, iteration, memory):
    # Update seed when called
    SEED + 1

    # Choose epsilon based on the iteration
    epsilon = get_epsilon_for_iteration(iteration)

    # Choose the action 
    if random.random() < epsilon:
        action = env.action_space.sample()
    else:
        action = choose_best_action(model, state)

    # Play one game iteration (note: according to the next paper, you should actually play 4 times here)
    [new_frame, reward, is_done, truncated, info] = env.step(action) # reward is based on color or brick destroyed
    memory.append([state, action, new_frame, reward, is_done])

    # Sample and fit, training model based on action for next iteration
    batch = memory.sample(1)

    # Encoding actions
    encodethis = np.array([item[1] for item in batch])
    encoded = np.empty((0,4), int)
    for item in encodethis:
      dummy = item
      item = []
      for i in range(4):
          if (dummy == i):
            item.append(1)
          else:
            item.append(0)
      encoded = np.append(encoded, np.array([item]), axis=0)

    start_state_list = np.array([item[0] for item in batch])
    next_states_list = np.array([item[2] for item in batch])
    rewards_list = np.array([item[3] for item in batch])
    terminal_list = np.array([item[4] for item in batch])
    fit_batch(model=model, gamma=0.99, start_states=start_state_list, actions=encoded, rewards=rewards_list, next_states=next_states_list, is_terminal=terminal_list)

    # returning necessary values for current iteration
    action_taken = action
    updated_state = new_frame
    terminal = is_done
    old_state = state
    return [updated_state, terminal]

# Main function definition
def main():

    # Notes on hyperparameters throughout the script, or can change to initialize them here. The uncommented values are to be implemented, except for n_episodes. 
    #minibatch_size = 32
    replay_memory = 1000000
    #agen_history_len = 4
    #discount = 0.99
    #repeat_frequency = 4
    #update_frequency = 4 # play 4 steps in the game for each Q-learning iteration
    #learning_rate = 0.00025 
    #gradient_momentum = 0.95
    #squared_gradient_momentum = 0.95
    #min_squared_gradient = 0.01
    #initial_exploration = 1
    #final_exploration = 0.1
    noopmax = 30 # do nothing actions to be performed at the start of an epsisode
    replay_start_zize = 50000 # uniform random policty is run for this number of frams before learning starts, and resulting experience is used to populate memory.
    #final_exploration_frame = 1000000
    n_episodes = 10 # change 100 - 1000 for better optimization
    show_every = 5

    # Initalize agent 
    agent = atari_model(N_ACTIONS)
    # Initialize memory class
    memory = RingBuf(replay_memory)
    # initialize random start state
    env = gym.make("BreakoutDeterministic-v4", render_mode="human")
    frame = env.reset()[0]

    # Key animation iteration
    for episode in range(n_episodes):
        print(f"Episode: {episode} Run")

        # Optional Preprocessing step here
        
        # Render
        env.render() 

        is_done = False
        iteration = 0
        while not is_done: # Perform a random action, returns the new frame, reward and whether the game is over
            # replace random action sample with q_iteration
            [frame, terminal] = q_iteration(env, agent, frame, iteration, memory) # setting new_state = to state to iterate
            # Render
            env.render()
            # update iteration value
            iteration += 1 
            # Look to end the loop
            is_done = terminal
        
        # Reset it, returns the starting frame
        frame = env.reset()[0]
        # Reset memory
        memory = RingBuf(replay_memory)
        
        env.close()

    agent.save('AtariModel')

# Call to main function
if __name__ == '__main__':
    main()
