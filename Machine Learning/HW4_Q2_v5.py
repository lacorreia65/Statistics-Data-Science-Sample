# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 23:51:09 2021

@author: LuisAlvaro
"""

## This is Q2.3 ****

import numpy as np
from qlearning import qlearn
from maze import MazeEnv, ProbabilisticMazeEnv
from plotting_utils import plot_steps_vs_iters, plot_several_steps_vs_iters, plot_policy_from_q

def main():
    
    num_iters = 200
    alpha = 1.0
    gamma = 0.9
    epsilon = 0.1
    max_steps = 100
    use_softmax_policy = True

    # Set the environment probability of random
    env_p_rand_list = [0.05, 0.1, 0.25, 0.5]

    steps_vs_iters_list = []
    for env_p_rand in env_p_rand_list:
        # Instantiate with ProbabilisticMazeEnv
        env = ProbabilisticMazeEnv()
        env.p_rand = env_p_rand

        # Note: We will repeat for several runs of the algorithm to make the result less noisy
        avg_steps_vs_iters = np.zeros(num_iters)
        for i in range(10):
            q_hat, steps_vs_iters = qlearn(env, num_iters, alpha, gamma, epsilon, max_steps, 
                               use_softmax_policy, init_beta=1.0, k_exp_sched=0.0, SSeed=579)
            avg_steps_vs_iters += steps_vs_iters
        avg_steps_vs_iters /= 10
        steps_vs_iters_list.append(avg_steps_vs_iters)

    # TODO: Plot the steps vs iterations
    label_list = ["env_random={}".format(env_p_rand) for env_p_rand in env_p_rand_list]
    plot_several_steps_vs_iters(steps_vs_iters_list, label_list, block_size=10)

if __name__ == "__main__":
    main()
