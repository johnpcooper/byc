import roi_management as rm

def run():
	
	cs = rm.Cell_Stacks()

	cropped_dsred_stacks, cropped_yfp_stacks, cropped_bf_stacks = cs.set_cropped_cells_lists(cs.cells_dfs)
	resized_dsred_stacks, resized_yfp_stacks, resized_bf_stacks = cs.set_resized_cells_lists(cropped_dsred_stacks, cropped_yfp_stacks, cropped_bf_stacks, cs.cells_dfs)
	cs.save_resized_stacks(resized_dsred_stacks, resized_yfp_stacks, resized_bf_stacks, cs.cells_dfs)
	otsu_thresh_dsred_stacks = cs.set_otsu_thresholded_cells_lists(resized_dsred_stacks, cs.cells_dfs)

if __name__ == "__main__":
	run()


