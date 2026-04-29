CC = g++
TARGET = a.out
SRCDIR = src
INCDIR = include
OBJDIR = obj
RESDIR = result
VTKDIR = vtk
GNUDIR = gnuplot
SRCS = $(wildcard $(SRCDIR)/*.cpp)
OBJS = $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SRCS))

all: $(TARGET)
	./$(TARGET)

$(TARGET): $(OBJS)
	$(CC) $^ -o $@

$(OBJS): $(OBJDIR)/%.o: $(SRCDIR)/%.cpp | $(OBJDIR)
	$(CC) -I$(INCDIR) -c $< -o $@

clean:
	rm -f $(OBJDIR)/*.o $(TARGET)

clean-all:
	rm -f $(OBJDIR)/*.o $(TARGET) $(RESDIR)/$(VTKDIR)/*.vtk