import .config
import .utils
import os

class Benchmark:

    def __init__(self, test_name, data=None):
        self.test_name = test_name
        self.json_file = os.path.join(config.benchmark_location, self.test_name + '.json')
        if data is None:
            self.data = self.read_json(self.test_name)
        else:
            self.data = data
            
    def reset():
        with open(self.json_file, 'w') as outfile:
            json.dump(self.data, outfile)        
            
    def get():
        with open('data.txt') as infile:
            data = json.load(infile)
        self.data = data
