
# Subcarrier Assignment, User Matching, and Power Allocation for Weighted Sum-Rate Maximization with RSMA

This repository contains simulations and analysis related to the proposed methodology in the paper:

**"Subcarrier Assignment, User Matching, and Power Allocation for Weighted Sum-Rate Maximization with RSMA"**

The code implements key algorithms, presents numerical results, and validates the theoretical findings discussed in the paper. For detailed explanations and theoretical foundations, please refer to the original article.

### Citation:
DE CASTRO, Jonas Alves; LIMA, Francisco Rafael Marques. *Subcarrier assignment, user matching and power allocation for weighted sum-rate maximization with RSMA*. In: 2022 Global Information Infrastructure and Networking Symposium (GIIS). IEEE, 2022. p. 20-24.

---

## Table of Contents

- [Introduction](#introduction)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Results](#results)
- [Acknowledgments](#acknowledgments)

---

## Introduction

This repository provides the code and tools to replicate the simulations and analyze the methodology described in the paper **"Subcarrier Assignment, User Matching, and Power Allocation for Weighted Sum-Rate Maximization with RSMA"**. The focus is on the following key aspects:

1. **Subcarrier Assignment**: Efficiently assigning subcarriers to users for optimal performance.
2. **User Matching**: Selecting pairs of users to maximize overall system performance.
3. **Power Allocation**: Allocating power across users and subcarriers to maximize the weighted sum-rate.

The code validates the theoretical results and provides numerical comparisons to demonstrate the efficacy of the proposed approach.

---

## Requirements

To run the simulations, the following software and toolboxes are required:

- **MATLAB** (or compatible environment)
- **Optimization Toolbox**: CPLEX is Required for solving optimization problems.
- **Communications System Toolbox**: Required for simulating communication systems.

---

## Installation

1. **Clone the repository**:

    ```bash
    git clone https://github.com/yourusername/repository.git
    ```

2. **Navigate to the repository directory**:

    ```bash
    cd repository
    ```

3. **Ensure that the required toolboxes are installed** in your MATLAB environment.

---

## Usage

1. **Run the main simulation script**:

    Navigate to the folder where the repository is located and execute the main simulation file, `main_simulation.m`:

    ```matlab
    main.m
    ```

2. **Configuration**: 
    - Before running, you may need to configure simulation parameters (e.g., number of users, subcarriers, channel conditions, etc.) in the `config.m` file.
    - The simulation script will generate results such as weighted sum-rate values and performance comparisons.

3. **Results**:
    - Results will be saved in the `results/` directory, which contains numerical output and figures.
    - You can analyze the performance and compare it with theoretical predictions discussed in the paper.

---

## Results

The simulation provides various performance metrics, including:

- Weighted sum-rate maximization performance.
- Comparison between the proposed methodology and baseline methods.
- Visualization of results in graphs and tables.

---

## Acknowledgments

- The methodology and algorithms in this repository are based on the work of **Jonas Alves de Castro** and **Francisco Rafael Marques Lima**.
- This project was developed as part of research efforts to validate and implement the results presented in their paper.

---
