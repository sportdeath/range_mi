class GridExplorer {

  GridExplorer(const GridMapper * mapper);

  /**
   * Each of these methods returns a path
   * that will .
   */
  std::vector<unsigned int> path_to_nearest_frontier();
  std::vector<unsigned int> path_to_nearest_of_k_max_mi(unsigned int k, map * );
  std::vector<unsigned int> path_maximizing_mi(bool per_length);

  private:
    const GridMapper * mapper;
}
