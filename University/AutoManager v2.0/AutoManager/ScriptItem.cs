using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace AutoManager {
	class ScriptItem {
		private string program, typeCommand, command;

		public ScriptItem(string program,string typeCommand, string command) {
			this.program = program;
			this.typeCommand = typeCommand;
			this.command = command;
		}

		public ScriptItem() {
			this.program = "";
			this.typeCommand = "";
			this.command = "";
		}

		public string Program {
			get { return program; }
			set { program = value; }
		}

		public string TypeCommand {
			get { return typeCommand; }
			set { typeCommand = value; }
		}

		public string Command {
			get { return command; }
			set { command = value; }
		}
	}
}
