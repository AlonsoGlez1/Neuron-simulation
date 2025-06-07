def calculate_mean_average_in_range(file_path, start_row, end_row):
    """
    Calculate the average of the 'mean' column for rows between start_row and end_row
    in a tab-separated file with specific headers.
    """
    try:
        with open(file_path, 'r') as file:
            lines = [line.strip() for line in file.readlines() if line.strip()]
            
            # Find header line
            header_line = None
            header_line_num = -1
            for i, line in enumerate(lines):
                if line.startswith('timeSteps\tmean\tstd_dev\tsem'):
                    header_line = line
                    header_line_num = i
                    break
            
            if not header_line:
                raise ValueError("Header line 'timeSteps\tmean\tstd_dev\tsem' not found")
            
            # Get column indices
            headers = header_line.split('\t')
            try:
                time_index = headers.index('timeSteps')
                mean_index = headers.index('mean')
            except ValueError as e:
                raise ValueError(f"Required column not found: {str(e)}")
            
            # Process data rows
            data_rows = lines[header_line_num+1:]
            total_rows = len(data_rows)
            
            # Validate row range
            if start_row < 1 or end_row > total_rows:
                raise ValueError(f"Row range must be between 1 and {total_rows}")
            if start_row > end_row:
                raise ValueError("Start row must be ≤ end row")
            
            # Extract values in range
            selected_values = []
            for i in range(start_row-1, end_row):
                if i >= len(data_rows):
                    break
                
                parts = data_rows[i].split('\t')
                if len(parts) <= mean_index:
                    continue
                
                try:
                    time_step = int(parts[time_index])
                    mean_value = float(parts[mean_index])
                    selected_values.append((time_step, mean_value))
                except ValueError:
                    continue
            
            if not selected_values:
                raise ValueError("No valid numeric data found in selected range")
            
            # Calculate average
            mean_values = [val for _, val in selected_values]
            average = sum(mean_values) / len(mean_values)
            
            return {
                'average': average,
                'time_range': (selected_values[0][0], selected_values[-1][0]),
                'values_used': len(selected_values),
                'all_values': selected_values
            }
            
    except FileNotFoundError:
        raise FileNotFoundError(f"File not found: {file_path}")
    except Exception as e:
        raise Exception(f"An error occurred: {str(e)}")

def main():
    print("Mean Value Range Calculator")
    print("--------------------------")
    print("For tab-separated files with header: timeSteps\tmean\tstd_dev\tsem\n")
    
    file_path = input("Enter file path: ").strip()
    
    try:
        # First scan to determine available rows
        with open(file_path, 'r') as f:
            lines = [line.strip() for line in f if line.strip()]
        
        # Find header
        header_pos = -1
        for i, line in enumerate(lines):
            if line.startswith('timeSteps\tmean\tstd_dev\tsem'):
                header_pos = i
                break
        
        if header_pos == -1:
            print("\nERROR: File must contain header: timeSteps\tmean\tstd_dev\tsem")
            return
        
        data_rows = lines[header_pos+1:]
        if not data_rows:
            print("\nERROR: No data rows found after header")
            return
        
        # Show basic info
        print(f"\nFound {len(data_rows)} data rows")
        print("First row:", data_rows[0])
        if len(data_rows) > 1:
            print("Last row:", data_rows[-1])
        
        # Get row range
        while True:
            try:
                start_row = int(input(f"\nEnter start row (1-{len(data_rows)}): "))
                end_row = int(input(f"Enter end row ({start_row}-{len(data_rows)}): "))
                
                if 1 <= start_row <= end_row <= len(data_rows):
                    break
                print("Invalid range. Try again.")
            except ValueError:
                print("Please enter valid numbers.")
        
        # Calculate average
        result = calculate_mean_average_in_range(file_path, start_row, end_row)
        
        # Display results
        print("\nRESULTS")
        print("-------")
        print(f"Time steps: {result['time_range'][0]} to {result['time_range'][1]}")
        print(f"Rows averaged: {result['values_used']}")
        print(f"Mean value average: {result['average']:.6f}")
        
        # Show sample values
        print("\nFirst 3 values:")
        for time, val in result['all_values'][:3]:
            print(f"Time {time}: {val:.6f}")
        
        if len(result['all_values']) > 3:
            print("\nLast 3 values:")
            for time, val in result['all_values'][-3:]:
                print(f"Time {time}: {val:.6f}")
        
        print("\nDone.")
        
    except Exception as e:
        print(f"\nERROR: {str(e)}")

if __name__ == "__main__":
    main()