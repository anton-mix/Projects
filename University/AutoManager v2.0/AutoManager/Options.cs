using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace AutoManager {
	class Options {
		List<string> childWindows, ellsements;
		string counterText;
		int counterStart, counterStep;

		public Options() {
			childWindows = new List<string>();
			ellsements = new List<string>();
			counterText = "";
			counterStart = counterStep = 1;
		}

		public List<string> ChildWindows {
			get { return childWindows; }
			set { childWindows = value; }
		}

		public List<string> Ellements {
			get { return ellsements; }
			set { ellsements = value; }
		}

		public string CounterText {
			get { return counterText; }
			set { counterText = value; }
		}

		public int CounterStart {
			get { return counterStart; }
			set { counterStart = value; }
		}

		public int CounterStep {
			get { return counterStep; }
			set { counterStep = value; }
		}

	}
}
