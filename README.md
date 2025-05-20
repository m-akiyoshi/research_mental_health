# Mental Health Therapy Recommendation System

## Overview

This AI-powered system provides personalized, evidence-based therapy recommendations based on an individual's mental health profile. By analyzing symptoms, history, traits, and goals, the system suggests appropriate therapies and explains its recommendations with supporting clinical literature.

## System Architecture

The system follows a multi-step pipeline:

1. **Data Collection**

   - Gathers research papers, therapy guidelines, and case studies
   - Focuses on evidence-based mental health interventions

2. **Text Processing**

   - Cleans and prepares text data
   - Implements NLP pipeline for high-quality input

3. **Knowledge Extraction**

   - Identifies key entities (symptoms, therapies, outcomes)
   - Extracts relationships between entities
   - Uses REBEL model for relation extraction

4. **Knowledge Graph**

   - Constructs graph representation of mental health concepts
   - Maps relationships between symptoms, therapies, and outcomes

5. **Pattern Recognition**

   - Generates graph embeddings
   - Identifies patterns and similarities
   - Clusters related concepts

6. **User Profiling**

   - Processes user input
   - Creates query graph from user profile

7. **Recommendation Engine**

   - Matches user profile with therapy options
   - Uses graph traversal for personalized recommendations

8. **Explanation Generation**
   - Provides rationale for recommendations
   - Cites supporting research
   - Generates interpretable output

## Setup and Installation

### Prerequisites

- Python 3.10 or higher
- pip (Python package manager)

### Installation

1. Clone the repository:

   ```bash
   git clone [repository-url]
   cd research_mental_health
   ```

2. Create and activate a virtual environment:

   ```bash
   python -m venv venv310
   source venv310/bin/activate  # On Windows: venv310\Scripts\activate
   ```

3. Install required packages:
   ```bash
   pip install -r requirements.txt
   ```

## License

This project is licensed under the MIT License - see the LICENSE file for details.
