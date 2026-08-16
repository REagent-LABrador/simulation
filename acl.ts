import type { ACL } from "@/lib/access.ts";

// Who may call this agent through the router. Enforcement starts once auth
// is wired into the eve router (see lib/access.ts).
//
// Public on purpose: nothing populates `ctx.session.auth` yet, and an
// unresolved caller fails closed — so a restricted ACL here would hide this
// tool from *every* router caller (HTTP, Slack, all of them) while
// `bun run console` kept working, which makes the gap easy to miss.
// Restrict once /managed-agent-setup has wired auth.
export const acl: ACL = { public: true };
// Restricted instead:
// export const acl: ACL = { principals: ["org_acme"] };
