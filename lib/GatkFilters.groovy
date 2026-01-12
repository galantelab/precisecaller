class GatkFilters {

	// Parse the input (String or List) into a structured List of Maps
	public static List parse(Object input) {
		// Case 1: User passed a List (from YAML/JSON params-file)
		if (input instanceof List) {
			return input.findAll { it instanceof Map && it.name && it.expr }
		}

		// Case 2: User passed a String (from CLI)
		if (input instanceof String && input.trim()) {
			return input.split(/(?<!\|)\|(?!\|)/)
				.collect { entry ->
					def parts = entry.split(':', 2)
					def name = parts.size() > 0 ? parts[0].trim() : ""
					def expr = parts.size() > 1 ? parts[1].trim() : ""
					[ name: name, expr: expr ]
				}
				.findAll { it.name && it.expr }
		}

		return []
	}

	// Transform the List of Maps into GATK-compatible CLI arguments
	public static String stringify(List filters) {
		if (!filters) return ""

		// Using single quotes around values to protect JEXL symbols from the shell
		return filters.collect {
			"--filter-name '${it.name}' --filter-expression '${it.expr}'"
		}.join(' ')
	}

}
