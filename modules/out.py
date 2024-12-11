import pyfiglet

# Output file writer

class Output_Writer:
    def __init__(self, filename):
        self.filename = filename

    def write_header(self,header):
        """Write a header to output file"""
        with open(self.filename, "w") as file:
             file.write(header + "\n")
    
    def write_subheader(self,subheader):
        """Writes a subheader for a new section"""
        with open(self.filename, "a") as file:
            file.write(subheader + "\n")
    
    def write_data(self,data):
        with open(self.filename, "a") as file:
            for line in data:
                file.write(line + "\n")

    def clear_file(self):
        """Clears the output of a current file"""
        open(self.filename, "w").close()

    def write_line(self,line):
        with open(self.filename, "a") as file:
            file.write(line + "\n")

    header = pyfiglet.figlet_format("GRAPH-IC", font="slant",width=110,justify="Center")
    


