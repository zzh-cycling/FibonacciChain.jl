JL = julia --project

default: init update test

init:
	$(JL) -e 'using Pkg; Pkg.precompile();'

init-docs:
	julia --project=docs -e 'using Pkg; Pkg.precompile();'

update:
	$(JL) -e 'using Pkg; Pkg.update(); Pkg.precompile();'

test:
	$(JL) -e 'using Pkg; Pkg.test()'

coverage:
	$(JL) -e 'using Pkg; Pkg.test(; coverage=true)'

serve:
	$(JL) -e 'using Pkg; Pkg.activate("docs"); using LiveServer; servedocs(;skip_dirs=["docs/src/assets", "docs/src/generated"])'

clean:
	rm -rf docs/build
	find . -name "*.cov" -type f -print0 | xargs -0 /bin/rm -f

# Get version from Project.toml
VERSION := $(shell grep '^version' Project.toml | sed 's/version = "\(.*\)"/\1/')

# Release workflow: tag, build docs, and register
release:
	@echo "Creating git tag v$(VERSION)..."
	git add -A
	git commit -m "Release v$(VERSION)" || true
	git tag -a "v$(VERSION)" -m "Release v$(VERSION)"
	git push origin dev
	git push origin "v$(VERSION)"
	@echo "Release v$(VERSION) complete. Register at JuliaRegistries if needed."

docs:
	julia --project=docs docs/make.jl

publish:
	@echo "To register this package, go to:"
	@echo "  https://github.com/JuliaRegistries/General"
	@echo "  or use @JuliaRegistrator bot: comment '@JuliaRegistrator register' on a commit/release"
	@echo ""
	@echo "Current version: v$(VERSION)"

bump-patch:
	@echo "Bumping patch version..."
	$(JL) -e 'using Pkg; lines = readlines("Project.toml"); \
		for (i, line) in enumerate(lines); \
			if startswith(line, "version"); \
				v = match(r"\"(\d+)\.(\d+)\.(\d+)\"", line); \
				newv = "$$(parse(Int, v[1])).$$(parse(Int, v[2])).$$(parse(Int, v[3])+1)"; \
				lines[i] = "version = \"$$newv\""; \
				println("New version: $$newv"); \
			end; \
		end; \
		write("Project.toml", join(lines, "\n") * "\n")'

bump-minor:
	@echo "Bumping minor version..."
	$(JL) -e 'using Pkg; lines = readlines("Project.toml"); \
		for (i, line) in enumerate(lines); \
			if startswith(line, "version"); \
				v = match(r"\"(\d+)\.(\d+)\.(\d+)\"", line); \
				newv = "$$(parse(Int, v[1])).$$(parse(Int, v[2])+1).0"; \
				lines[i] = "version = \"$$newv\""; \
				println("New version: $$newv"); \
			end; \
		end; \
		write("Project.toml", join(lines, "\n") * "\n")'

bump-major:
	@echo "Bumping major version..."
	$(JL) -e 'using Pkg; lines = readlines("Project.toml"); \
		for (i, line) in enumerate(lines); \
			if startswith(line, "version"); \
				v = match(r"\"(\d+)\.(\d+)\.(\d+)\"", line); \
				newv = "$$(parse(Int, v[1])+1).0.0"; \
				lines[i] = "version = \"$$newv\""; \
				println("New version: $$newv"); \
			end; \
		end; \
		write("Project.toml", join(lines, "\n") * "\n")'

.PHONY: init test coverage serve clean update release docs publish bump-patch bump-minor bump-major
