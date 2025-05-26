import numpy as np
from scipy.integrate import simpson
from scipy.spatial.distance import euclidean
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#__after imporitng Libs 
class MobiusStrip:
    def __init__(self, radius, width, resolution):
      self.radius = radius  # Main radius of the Möbius strip
      self.width = width    # Width of the strip
      self.resolution = resolution  # Number of points for discretization
      self.u = np.linspace(0, 2 * np.pi, resolution)  # Parameter for the loop around the strip
      self.v = np.linspace(-width / 2, width / 2, resolution)  # Parameter across the width
      self.u_grid, self.v_grid = np.meshgrid(self.u, self.v)  # Create 2D grid of parameters
      self.x, self.y, self.z = self.compute_mesh()  # Compute 3D coordinates

    def compute_mesh(self):
      u = self.u_grid  # Get the u parameter grid
      v = self.v_grid  # Get the v parameter grid
      R = self.radius  # Main radius
      
      # Parametric equations for Mobius strip:
      x = (R + v * np.cos(u / 2)) * np.cos(u)  # x-coordinate
      y = (R + v * np.cos(u / 2)) * np.sin(u)  # y-coordinate
      z = v * np.sin(u / 2)  # z-coordinate
      
      return x, y, z

    def surface_area(self):
      # Compute partial derivatives of each coordinate with respect to u and v
      dx_du = np.gradient(self.x, axis=1)  # dx/du
      dx_dv = np.gradient(self.x, axis=0)  # dx/dv
      dy_du = np.gradient(self.y, axis=1)  # dy/du
      dy_dv = np.gradient(self.y, axis=0)  # dy/dv
      dz_du = np.gradient(self.z, axis=1)  # dz/du
      dz_dv = np.gradient(self.z, axis=0)  # dz/dv
      
      # Compute cross product of tangent vectors (dr/du × dr/dv)
      cross_x = dy_du * dz_dv - dz_du * dy_dv  # x-component of cross product
      cross_y = dz_du * dx_dv - dx_du * dz_dv  # y-component
      cross_z = dx_du * dy_dv - dy_du * dx_dv  # z-component
      
      # Magnitude of the cross product gives differential area element
      dA = np.sqrt(cross_x**2 + cross_y**2 + cross_z**2)
      
      # Double integral using Simpson's rule to compute total area
      area = simpson([simpson(row, self.u) for row in dA], self.v)
      return area

    def edge_length(self):
      # Extract top and bottom edges (actually same edge for Möbius strip)
      top_edge = np.array([self.x[-1], self.y[-1], self.z[-1]]).T
      bottom_edge = np.array([self.x[0], self.y[0], self.z[0]]).T
      
      def compute_length(edge):
          # Sum Euclidean distances between consecutive points
          return sum(
              euclidean(edge[i], edge[i + 1])
              for i in range(len(edge) - 1)
          )
      
      # Return total length (though for Möbius strip it's actually one edge)
      return compute_length(top_edge) + compute_length(bottom_edge)

    def plot(self):
        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(self.x, self.y, self.z, rstride=1, cstride=1,
                        cmap='viridis', edgecolor='k', linewidth=0.1)
        ax.set_title("Möbius Strip")
        plt.tight_layout()
        plt.show()
#__class and methods are created

if __name__ == "__main__":
    strip = MobiusStrip(radius=5, width=1, resolution=200)
    
    print("Computing properties...")
    print(f"Surface Area ≈ {strip.surface_area():.4f} units²")
    print(f"Edge Length ≈ {strip.edge_length():.4f} units")
    
    strip.plot()



