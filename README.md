# Dengue Transmission and Human Mobility Research Project

## Project Overview

This project is a comprehensive study on the relationship between dengue transmission dynamics and human mobility, employing mathematical modeling, machine learning, and spatial analysis methods to investigate dengue transmission patterns and risk assessment across different regions. The project primarily focuses on dengue outbreak data from the Xishuangbanna region, integrating climate factors and human mobility data for modeling analysis.

## Project Structure

### Data Processing Module

#### `01_get_dengue_data.R`
- **Function**: Preprocessing and cleaning of dengue case data
- **Main Operations**:
  - Read 2023 dengue case data
  - Data format conversion and date processing
  - Aggregate case counts by date
  - Generate complete time series data (fill missing dates)
- **Output**: Cleaned dengue case time series data

### SEI-SEIR Model Module

#### `02_1_SEI-SEIR_simulation_setup.R`
- **Function**: Basic parameter setup for SEI-SEIR model
- **Main Content**:
  - Load climate data
  - Set study location (Xishuangbanna)
  - Define population parameters (population size, birth rate, death rate)
  - Set model time step and migration parameters

#### `02_2_SEI-SEIR_model_THR.R`
- **Function**: Core implementation of Temperature-Humidity-Rainfall dependent SEI-SEIR model
- **Main Content**:
  - Define SEI-SEIR differential equation system
  - Implement temperature-dependent entomological parameter functions
  - Include Briere function, quadratic function and other biological response functions
  - Model considers three mosquito developmental stages and four human health states

#### `02_3_SEI-SEIR_simulate_THR_model.R`
- **Function**: Run SEI-SEIR model simulations
- **Main Operations**:
  - Perform model simulations using different rainfall functions
  - Compare three rainfall response modes: Briere, quadratic, and inverse functions
  - Generate and save model simulation results

#### `02_4_SEI-SEIR_parameter estimation-gibbs.R`
- **Function**: Parameter estimation using Gibbs sampling
- **Main Content**:
  - Implement Bayesian parameter estimation framework
  - Define loss function and likelihood function
  - Estimate model parameters using MCMC methods
  - Parameter uncertainty quantification

#### `02_6_SEIR and SEI-SEIR_parameter estimation.R`
- **Function**: Parameter estimation comparison between SEIR and SEI-SEIR models

### Standard SEIR Model Module

#### `03_1_SEIR_model_simulation.R`
- **Function**: Implementation and simulation of standard SEIR model
- **Main Content**:
  - Implement classical SEIR epidemic model
  - Include seasonal transmission rate variations
  - Model parameter calibration and simulation execution

#### `03_IR(t) calculation.R`
- **Function**: Calculate infection rate time series IR(t)
- **Main Operations**:
  - Process model simulation results
  - Calculate daily new infections
  - Analyze transmission risk combined with traveler data
  - Data visualization and result presentation

### Human Mobility Analysis Module

#### `04_01_migration matrix and simulation.R`
- **Function**: Population migration matrix construction and simulation
- **Main Content**:
  - Process Baidu migration data
  - Construct emigration index and migration proportion matrix
  - Analyze population flow patterns between Xishuangbanna and other cities
  - Population flow simulation and visualization

### Machine Learning Module

#### `05_01_ML data collecting.R`
- **Function**: Machine learning data collection and organization
- **Main Operations**:
  - Integrate traveler case migration simulation results
  - Calculate average migration case numbers and first arrival times
  - Merge local and travel transmission data
  - Prepare machine learning feature data

#### `05_02_ML EDA.R`
- **Function**: Machine learning exploratory data analysis
- **Main Content**:
  - Data distribution visualization
  - Variable correlation analysis
  - Data preprocessing and standardization
  - Feature engineering and variable selection

#### `05_03_validation.R`
- **Function**: Model validation and result comparison
- **Main Operations**:
  - Compare model predictions with actual data
  - Calculate correlation coefficients and statistical significance
  - Provincial data aggregation and validation
  - Result visualization

### Subsequent Risk Assessment Module

#### `06_01_subsequent risk data collection.R`
- **Function**: Subsequent risk assessment data collection
- **Main Content**:
  - Collect national climate data (temperature, rainfall)
  - Data format conversion and cleaning
  - Prepare input data for risk assessment models

#### `06_02_SEI-SEIR_model_THR with specific imported cases.R`
- **Function**: SEI-SEIR model considering specific imported cases

#### `06_02_subsequent risk setup and simulation.R`
- **Function**: Subsequent risk assessment setup and simulation

## Workflow

### Phase 1: Data Preparation
1. **Data Collection** (`01_get_dengue_data.R`)
   - Dengue case data preprocessing
   - Time series data construction

### Phase 2: Model Development
2. **SEI-SEIR Model Construction** (`02_*` series files)
   - Model parameter setup
   - Core model implementation
   - Model simulation execution
   - Parameter estimation and calibration

3. **SEIR Model Comparison** (`03_*` series files)
   - Standard SEIR model implementation
   - Infection rate calculation
   - Model result comparison

### Phase 3: Spatial Analysis
4. **Human Mobility Analysis** (`04_01_migration matrix and simulation.R`)
   - Migration data processing
   - Mobility matrix construction
   - Spatial transmission simulation

### Phase 4: Machine Learning
5. **Predictive Model Development** (`05_*` series files)
   - Feature engineering
   - Model training
   - Result validation

### Phase 5: Risk Assessment
6. **Subsequent Risk Assessment** (`06_*` series files)
   - Risk data collection
   - Risk model simulation
   - Prediction and assessment

## Main Dependencies

```r
# Data processing
library(dplyr)
library(reshape2)
library(lubridate)
library(readxl)

# Numerical computation and modeling
library(deSolve)      # Differential equation solving
library(mcmc)         # MCMC sampling
library(truncnorm)    # Truncated normal distribution

# Machine learning
library(caret)        # Model training and validation
library(corrplot)     # Correlation visualization
library(car)          # Variance inflation factor

# Spatial analysis
library(sf)           # Spatial data processing
library(sp)           # Spatial data structures
library(spdep)        # Spatial dependence analysis
library(rgdal)        # Geospatial data abstraction library

# Visualization
library(ggplot2)      # Data visualization
library(cowplot)      # Plot composition
library(ggspatial)    # Spatial data visualization
```

## Data Requirements

### Input Data
- **Dengue case data**: 2023 case onset dates and address information
- **Climate data**: Temperature, humidity, rainfall time series
- **Human mobility data**: Baidu migration index and inter-city flow proportions
- **Demographic data**: Regional population size, birth rates, death rates

### Output Data
- **Model simulation results**: Time series predictions from SEI-SEIR and SEIR models
- **Parameter estimation results**: Posterior distributions of model parameters
- **Risk assessment results**: Dengue transmission risk scores for different regions
- **Machine learning predictions**: Multi-factor based transmission risk predictions

## Usage Instructions

1. **Environment Setup**: Ensure all required R packages are installed
2. **Data Preparation**: Place raw data in the corresponding `data/` directory
3. **Sequential Execution**: Run scripts in order according to file numbering
4. **Result Review**: Check result files in the `output/` directory

## Important Notes

- Ensure correct working directory setup
- Some scripts require extended computation time (especially MCMC sampling)
- Some file paths may need adjustment based on actual conditions
- Recommend checking data file integrity before execution

## Research Significance

This project provides a comprehensive analytical framework for transmission prediction and risk assessment of dengue and other vector-borne diseases by integrating infectious disease dynamics modeling, spatial analysis, and machine learning methods. The research results can provide scientific evidence for public health decision-making and epidemic prevention and control.