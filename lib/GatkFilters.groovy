class GatkFilters {

	// Parse the input (String or Map)
	public static Map parse(Object input) {
		// Case 1: User passed a Map (from YAML/JSON params-file)
		if (input instanceof Map) {
			return input.findAll { k, v -> k && v }
		}

		// Case 2: User passed a String (from CLI)
		if (input instanceof String && input.trim()) {
			// Remove newlines/extra whitespace for multi-line config support
			String clean_input = input.replaceAll(/[\n\r\s]+/, " ")

			return clean_input.split(/(?<!\|)\|(?!\|)/)
				.collectEntries { entry ->
					def parts = entry.split(':', 2)
					if (parts.size() == 2) {
						return [parts[0].trim(), parts[1].trim()]
					}
					return [:]
				}
				.findAll { k, v -> k && v }
		}

		return [:]
	}

}
