
import os
import matplotlib.pyplot as plt

directory = 'Parker_EPR/ONME'
shortest_range = None
for filename in os.listdir(directory):
    if filename.endswith('.txt'):
        with open(os.path.join(directory, filename), 'r') as file:
            lines = file.readlines()
        data = [[float(value) for value in line.split()] for line in lines]
        x_range = max(row[0] for row in data) - min(row[0] for row in data)
        if shortest_range is None or x_range < shortest_range:
            shortest_range = x_range
# Create a dictionary to keep track of the number of plots for each complex
complex_counts = {}
loop = 0
for filename in os.listdir(directory):
    # Only process .txt files
    if filename.endswith('.txt'):
        # Open the file
        with open(os.path.join(directory, filename), 'r') as file:
            # Read the data
            lines = file.readlines()
        parts = filename.split('complex')
        if len(parts) > 1:
            label = parts[1].split('_')[0]
        else:
            label = 'No complex'
        # Increment the count for this complex
        if label in complex_counts:
            complex_counts[label] += 1
        else:
            complex_counts[label] = 1
        # Add the count to the label
        label += str(complex_counts[label])
        label += '='
        data = [[float(value) for value in line.split()] for line in lines]
        five_percent = int(len(data) * 0.05)
        data = data[five_percent:-five_percent]
        # Determine the middle range of x-values
        middle_x_start = len(data)// 4  # Start at the 1st quartile
        middle_x_end = 3*middle_x_start  # End at the 3rd quartile
        # Find the indices of the x-values that fall within this range
        middle_x_values = []
        middle_y_values = []
        i = middle_x_start
        while i < middle_x_end:
            middle_x_values.append(data[i][0])
            middle_y_values.append(data[i][1])
            i += 1
        max_y = max(middle_y_values)
        min_y = min(middle_y_values)
        meany = sum(middle_y_values) / len(middle_y_values)
        # Find the x-values corresponding to these y-values
        max_x = middle_x_values[middle_y_values.index(max_y)]
        min_x = middle_x_values[middle_y_values.index(min_y)]
        # Find the difference in x-values
        diff = min_x - max_x
        label +=str(round(diff,4))
        offset = loop * 600
        loop +=1
        if loop == 1:
            plt.plot([row[0] for row in data], [10*(row[1] - meany) + offset for row in data], label=label)
        elif loop == 2:
            plt.plot([row[0] for row in data], [2*(row[1] - meany) + offset for row in data], label=label)
        else:
            plt.plot([row[0] for row in data], [row[1] - meany + offset for row in data], label=label)
            
            
            
plt.xlim(12140,12260)
plt.xlabel('Magnetic Field Strength [Gauss]')
plt.ylabel('Intensity [a.u.]')
plt.yticks([])
plt.title('EPR Spectra')
plt.legend(title='Gauss (x) of ymin - ymax', loc='lower left')
plt.show()

exit()
directory = 'Parker_EPR'
shortest_range = None
for filename in os.listdir(directory):
    if filename.endswith('.txt'):
        with open(os.path.join(directory, filename), 'r') as file:
            lines = file.readlines()
        data = [[float(value) for value in line.split()] for line in lines]
        x_range = max(row[0] for row in data) - min(row[0] for row in data)
        if shortest_range is None or x_range < shortest_range:
            shortest_range = x_range
            
loop = 0
for filename in os.listdir(directory):
    # Only process .txt files
    if filename.endswith('.txt'):
        # Open the file
        with open(os.path.join(directory, filename), 'r') as file:
            # Read the data
            lines = file.readlines()
        parts = filename.split('complex')
        if len(parts) > 1:
            label = parts[1].split('_')[0]
        else:
            label = 'No complex'

        # Extract the desired part of the filename
        start = 36
        end = filename.find('_', start)  # Find the next underscore after the 36th character
        if end == -1:  # If there is no underscore after the 36th character, use the rest of the filename
            end = len(filename)
        plot_title = filename[start:end]

        data = []
        # ... rest of the code ...
        for line in lines:
            # Split the line into values
            values = line.split()
            # Convert the values to floats and add them to the data list
            data.append([float(value) for value in values])
        five_percent = int(len(data) * 0.05)
        data = data[five_percent:-five_percent]
        # Plot the data with an offset on the y-axis
        yoffset = loop * 4000  # Adjust the multiplier as needed
        loop +=1
        plt.plot([row[0] for row in data], [row[1] + yoffset for row in data], label=plot_title)

# ... rest of the code ...

# Add a legend and show the plot
plt.xlim(12000, 12500)
plt.legend()
plt.show()



exit()
# Directory containing the files
directory = 'Parker_EPR'
# Loop through all files in the directory
loop = 0
for filename in os.listdir(directory):
    loop += 1
    # Only process .txt files
    if filename.endswith('.txt'):
        # Open the file
        with open(os.path.join(directory, filename), 'r') as file:
            # Read the data
            lines = file.readlines()
        # Extract the desired part of the filename
        start = 36
        end = filename.find('_', start)  # Find the next underscore after the 36th character
        if end == -1:  # If there is no underscore after the 36th character, use the rest of the filename
            end = len(filename)
        plot_title = filename[start:end]

        data = []
        # ... rest of the code ...
        for line in lines:
            # Split the line into values
            values = line.split()
            # Convert the values to floats and add them to the data list
            data.append([float(value) for value in values])
        
        # Determine the middle range of x-values
        middle_x_start = len(data)// 4  # Start at the 1st quartile
        middle_x_end = 3*middle_x_start  # End at the 3rd quartile
        # Find the indices of the x-values that fall within this range
        middle_x_values = []
        middle_y_values = []
        i = middle_x_start
        while i < middle_x_end:
            middle_x_values.append(data[i][0])
            middle_y_values.append(data[i][1])
            i += 1
        max_y = max(middle_y_values)
        min_y = min(middle_y_values)
        
        # Find the x-values corresponding to these y-values
        max_x = middle_x_values[middle_y_values.index(max_y)]
        min_x = middle_x_values[middle_y_values.index(min_y)]
        # Find the difference in x-values
        diff = max_x - min_x
        print(filename)
        print(loop, diff)
        # Calculate the indices that correspond to the first 5% and last 5% of the list
        start_index = int(len(data) * 0.05)
        end_index = int(len(data) * 0.95)
        # Use list slicing to get the middle 90% of the list
        data = data[start_index:end_index]
        # Now you can use trimmed_data in place of data
        middle_x_values = []
        middle_y_values = []
        # Plot the data
        plt.figure()  # Create a new figure
        plt.plot([row[0] for row in data], [row[1] for row in data])
        plt.xlabel('Magnetic Field Strength [Gauss]')
        plt.ylabel('Intensity [a.u.]')
        plt.title(plot_title)  # Use the filename as the plot title
plt.show()  # Display all plots
exit()

# Open the file
with open('Parker_EPR/x224041104.DSC_x224041104_Gd_complexONMe_ESE_pi80_tau500_10K_9.534148 GHz.txt', 'r') as file:
    # Read the data
    lines = file.readlines()

data = []
for line in lines:
    # Split the line into values
    values = line.split()
    # Convert the values to floats and add them to the data list
    data.append([float(value) for value in values])

print(data[0])

#plot the data.
plt.plot(data[:][0], data[:][1])
plt.xlabel('B')
plt.ylabel('Frequency')
plt.title('EPR Spectrum')
plt.show()
