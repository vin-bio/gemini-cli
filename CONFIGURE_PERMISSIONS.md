# Configuring Gemini CLI for Maximum Action Permissions

This guide ensures that all actions the agentic system takes are always allowed without manual confirmations.

## 1. YOLO Mode (Recommended)

**YOLO mode** is the most comprehensive way to bypass all confirmations:

```bash
# Always run with YOLO mode
gemini --yolo
# or
gemini -y

# Create an alias for convenience
alias gemini-auto="gemini --yolo"
```

**In YOLO mode:**
- All tool confirmations are skipped
- File edits are automatically approved
- Shell commands execute without prompts
- Web fetches proceed automatically
- MCP tools run without user confirmation

## 2. Environment Variables

Set these in your `.bashrc` or `.env` file:

```bash
# API Authentication
export GEMINI_API_KEY="your-api-key-here"
export GEMINI_MODEL="gemini-2.5-pro"

# Disable sandboxing for maximum permissions
export GEMINI_SANDBOX=false

# For Vertex AI (alternative)
# export GOOGLE_API_KEY="your-google-api-key"
# export GOOGLE_GENAI_USE_VERTEXAI=true
# export GOOGLE_CLOUD_PROJECT="your-project-id"
```

## 3. Settings Configuration

Create/update `.gemini/settings.json`:

```json
{
  "comment": "Configuration for maximum permissions",
  
  "sandbox": false,
  "autoAccept": true,
  
  "mcpServers": {},
  "allowMCPServers": ["*"],
  "excludeMCPServers": [],
  
  "coreTools": [
    "list_directory",
    "read_file", 
    "write_file",
    "edit",
    "grep",
    "glob",
    "run_shell_command",
    "web_fetch",
    "web_search",
    "read_many_files",
    "memory"
  ],
  "excludeTools": [],
  
  "fileFiltering": {
    "respectGitIgnore": false,
    "enableRecursiveFileSearch": true
  }
}
```

## 4. Approval Modes

The system has three approval modes:

- **DEFAULT**: Requires confirmation for most actions
- **AUTO_EDIT**: Auto-approves file editing actions
- **YOLO**: Auto-approves ALL actions (recommended)

### Runtime Controls

During a session, you can toggle modes:
- `Ctrl+Y`: Toggle YOLO mode on/off
- `Shift+Tab`: Toggle AUTO_EDIT mode

## 5. Sandboxing Configuration

Disable sandboxing to allow maximum system access:

```bash
# Command line
gemini --no-sandbox

# Environment variable
export GEMINI_SANDBOX=false

# Settings file
{
  "sandbox": false
}
```

If you must use sandboxing, use the most permissive profile:

```bash
# macOS Seatbelt - most permissive
export SEATBELT_PROFILE=permissive-open
```

## 6. MCP Server Trust

For external MCP servers, configure them as trusted:

```json
{
  "mcpServers": {
    "your-server": {
      "command": "your-command",
      "trust": true
    }
  },
  "allowMCPServers": ["*"]
}
```

## 7. Tool-Specific Configurations

### Shell Commands
Commands are automatically whitelisted after first approval. In YOLO mode, all commands run immediately.

### File Operations
- `write_file`: Auto-approved in YOLO/AUTO_EDIT modes
- `edit`: Auto-approved in YOLO/AUTO_EDIT modes
- `read_file`: Generally no confirmation required

### Web Operations
- `web_fetch`: Auto-approved in YOLO/AUTO_EDIT modes
- `web_search`: Built-in tool, typically no confirmation

## 8. Complete Example Configuration

For a fully permissive setup:

**`.bashrc`:**
```bash
export GEMINI_API_KEY="your-key-here"
export GEMINI_MODEL="gemini-2.5-pro"
export GEMINI_SANDBOX=false
alias gemini="gemini --yolo"
```

**`.gemini/settings.json`:**
```json
{
  "sandbox": false,
  "autoAccept": true,
  "allowMCPServers": ["*"],
  "excludeTools": [],
  "fileFiltering": {
    "respectGitIgnore": false,
    "enableRecursiveFileSearch": true
  },
  "maxSessionTurns": 1000
}
```

## 9. Verification

To verify your configuration:

```bash
# Check if YOLO mode is active (look for red "YOLO mode" indicator)
gemini

# Test tool execution without confirmations
gemini -p "create a test file and write hello world to it"

# Verify sandbox status
gemini -p "run shell command: echo $SANDBOX"
```

## 10. Security Considerations

⚠️ **Important Security Notes:**

- YOLO mode disables all safety confirmations
- Shell commands will execute without review
- File modifications happen automatically
- Web requests are made without prompts
- Only use in trusted environments
- Consider using sandboxing in production environments
- Review the actions the AI takes regularly

## 11. Troubleshooting

If actions are still being blocked:

1. **Check approval mode**: Look for indicators in the CLI
2. **Verify environment variables**: `echo $GEMINI_SANDBOX`
3. **Review settings**: Check `.gemini/settings.json`
4. **Test incrementally**: Start with read-only operations
5. **Check logs**: Enable debug mode with `DEBUG=1 gemini`

## 12. Advanced Configuration

For even more control, you can:

- Disable telemetry: `"telemetry": {"enabled": false}`
- Extend session limits: `"maxSessionTurns": 10000`
- Configure tool summarization: `"summarizeToolOutput": {}`
- Set custom tool discovery: `"toolDiscoveryCommand": "custom-command"`

With these configurations, the Gemini CLI will operate with maximum autonomy and minimal user intervention. 