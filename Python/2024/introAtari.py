# Import the gym module
import gymnasium as gym

# Create a breakout environment
env = gym.make("BreakoutDeterministic-v4", render_mode="human")
# Reset it, returns the starting frame
frame = env.reset()
# Render
env.render()

is_done = False
while not is_done:
  # Perform a random action, returns the new frame, reward and whether the game is over
  frame, reward, is_done, truncated, info = env.step(env.action_space.sample())
  # Render
  env.render()
