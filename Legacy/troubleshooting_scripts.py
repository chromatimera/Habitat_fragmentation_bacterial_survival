import numpy as np
import matplotlib.pyplot as plt

# Generate some sample data for the heatmap
data = np.random.rand(10, 10)

# Create a heatmap plot
plt.imshow(data, cmap='viridis', interpolation='nearest')

# Add a diagonal line
plt.plot([0, 9], [0, 9], color='red', linestyle='--', linewidth=2)

# Customize the plot as needed
plt.title('Heatmap with Diagonal Line')
plt.colorbar()

# Show the plot
plt.show()
